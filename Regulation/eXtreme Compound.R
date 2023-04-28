

if(T){
  library(patchwork)
  library(ggalluvial)
  library(patchwork)
  library(reshape2)
  library(scales)
  library(ggpubr)
  library(pheatmap)
  library(data.table)
  library(reshape2)
  library(ggthemes)
  library(ggplot2)
  library(Seurat)
  library(cowplot)
  library(patchwork)
  library(Matrix)
  library(viridis)
  library(RColorBrewer)
  library(ggpubr)
  library(tictoc)
  library(RColorBrewer)
  library(corrplot)
  library(GeneOverlap)
  library(grid)
  library(gridExtra)
  library(igraph)
  library(ggrepel)
  library(tidyverse)
  library(paletteer)
  library(ggforce)
  library(ggrastr)
  library(rms)
  library(qusage)
  library(spiralize)
  library(doParallel)
  library(pROC)
  library(vegan)
  library(tidyverse)
  library(dplyr)
  library(data.table)
  library(circlize)
  library(ComplexHeatmap )
  library(WGCNA)
}



if (T) {
  dd1 <- list()
  for (i in list.files(path = "ref1/", pattern = "csv")) {
    dd1[[paste0("Cohort", gsub("GSE|_exprSet|\\.csv", "", i))]] <- fread(paste0("ref1/", i), data.table = F) %>% column_to_rownames("V1")
  }
  dd2 <- list()
  for (i in list.files(path = "ref2/", pattern = "csv")) {
    dd2[[paste0("Cohort", gsub("GSE|_targets|\\.csv", "", i))]] <- read.csv(paste0("ref2/", i), row.names = 1)
  }
}
APVS <- read.table("APVS gene.txt") %>% .[, 1]


dd1_APVSLevel <- lapply(dd1, function(x) {
  x <- x[ which(rownames(x) %in% gene), ]
  brca.dist <- vegdist(t(x),
                       method = "euclidean"
  )
  pc <- cmdscale(brca.dist,
                 eig = FALSE
  )
  geneprofile <- x %>% t()
  pre <- apply(
    geneprofile, 2,
    function(m) {
      df <- m * (pc[, 1] + pc[, 2])
    }
  )
  data <- data.frame(
    APVSLevel = apply(pre, 1, sum)
  )
})



####################################################
#######    Aberrant expression pattarn    ##########
####################################################

if (T) {
  enableWGCNAThreads()
  fpkm <- dd1[["Cohort 1"]]
  WGCNA_matrix <- t(fpkm[order(apply(fpkm, 1, mad), decreasing = T)[1:5000], ])
  datExpr0 <- WGCNA_matrix
  gsg <- goodSamplesGenes(datExpr0, verbose = 3)
  sampleTree <- hclust(dist(datExpr0), method = "average")
  clust <- cutreeStatic(sampleTree, cutHeight = 51, minSize = 10)
  keepSamples <- (clust == 1)
  datExpr <- datExpr0[keepSamples, ]
  datExpr <- datExpr0
  nGenes <- ncol(datExpr)
  nSamples <- nrow(datExpr)
  save(datExpr, nGenes, nSamples, file = "Step01-WGCNA_input.Rdata")
  powers <- c(c(1:10), seq(from = 12, to = 28, by = 2))
  sft <- pickSoftThreshold(datExpr, powerVector = powers, verbose = 5)

  softpower <- sft$powerEstimate
  ADJ <- abs(cor(datExpr, use = "p"))^softpower
  k <- as.vector(apply(ADJ, 2, sum, na.rm = T))

  adjacency <- adjacency(datExpr, power = softpower)
  TOM <- TOMsimilarity(adjacency)
  dissTOM <- 1 - TOM
  library(flashClust)
  geneTree <- flashClust(as.dist(dissTOM), method = "average")
  plot(geneTree, xlab = "", sub = "", main = "Gene clustering on TOM-based", labels = FALSE, hang = 0.04)

  minModuleSize <- 30
  dynamicMods <- cutreeDynamic(
    dendro = geneTree,
    distM = dissTOM,
    deepSplit = 2,
    pamRespectsDendro = FALSE,
    minClusterSize = minModuleSize
  )
  dynamicColors <- labels2colors(dynamicMods)
  plotDendroAndColors(
    dendro = geneTree,
    colors = dynamicColors,
    groupLabels = "Dynamic Tree Cut",
    dendroLabels = FALSE, hang = 0.03,
    addGuide = TRUE,
    main = "Gene dendrogram and module colors"
  )

  MEList <- moduleEigengenes(datExpr, colors = dynamicColors)
  MEs <- MEList$eigengenes
  MEDiss1 <- 1 - cor(MEs)
  METree1 <- flashClust(as.dist(MEDiss1), method = "average")
  MEDissThres <- 0.25
  merge <- mergeCloseModules(datExpr,
    dynamicColors,
    cutHeight = MEDissThres,
    verbose = 3
  )
  mergedColors <- merge$colors
  mergedMEs <- merge$newMEs

  moduleColors <- mergedColors
  colorOrder <- c("grey", standardColors(50))
  moduleLabels <- match(moduleColors, colorOrder) - 1
  MEs <- mergedMEs
  MEDiss2 <- 1 - cor(MEs)
  METree2 <- flashClust(as.dist(MEDiss2), method = "average")

  plotDendroAndColors(
    dendro = geneTree,
    colors = cbind(dynamicColors, mergedColors),
    groupLabels = c("Dynamic Tree Cut", "Merged Dynamics"),
    dendroLabels = FALSE,
    hang = 0.04,
    addGuide = TRUE,
    guideHang = 0.25,
    cex.colorLabels = 1,
    colorHeight = 0.2,
    colorHeightBase = 0.3,
    colorHeightMax = 0.6,
    main = "Gene Dendrogram and module colors"
  )

  write.table(table(moduleColors), "MEgeneCount.txt", quote = F, row.names = F)
  moduleColors <- mergedColors
  colorOrder <- c("grey", standardColors(50))
  moduleLabels <- match(moduleColors, colorOrder) - 1
  MEs <- mergedMEs
  MEs <- orderMEs(MEs)

  ## TOMplot
  softpower <- 10
  dissTOM <- 1 - TOMsimilarityFromExpr(datExpr, power = softpower)
  nSelect <- 500
  set.seed(10)
  select <- sample(nGenes, size = nSelect)
  selectTOM <- dissTOM[select, select]
  selectTree <- hclust(as.dist(selectTOM), method = "average")
  selectColors <- moduleColors[select]
  sizeGrWindow(9, 9)
  plotDiss <- selectTOM^7
  diag(plotDiss) <- NA
  TOMplot(plotDiss, selectTree, selectColors,
    main = "Network heatmap plot, selected genes"
  )
}


Clincal <- dd1_APVSLevel[["Cohort"]]
datTraits <- clinical
rownames(datTraits) <- rownames(clinical)
datTraits <- as.data.frame(do.call(cbind, lapply(datTraits, as.numeric)))
rownames(datTraits) <- rownames(clinical)
datTraits <- datTraits[row.names(datExpr), , drop = F]


sampleTree2 <- hclust(dist(datExpr), method = "average")
traitColors <- numbers2colors(datTraits, signed = FALSE)
plotDendroAndColors(sampleTree2,
  traitColors,
  groupLabels = names(datTraits),
  main = "Sample dendrogram and trait heatmap"
)


MEs <- orderMEs(MEs)
moduleTraitCor <- cor(MEs, datTraits, use = "p")
moduleTraitPvalue <- corPvalueStudent(moduleTraitCor, nSamples)

textMatrix <- paste(signif(moduleTraitCor, 2), "\n(", signif(moduleTraitPvalue, 1), ")", sep = "")
dim(textMatrix) <- dim(moduleTraitCor)
labeledHeatmap(
  Matrix = moduleTraitCor,
  xLabels = colnames(datTraits),
  yLabels = names(MEs),
  ySymbols = names(MEs),
  colorLabels = FALSE,
  colors = blueWhiteRed(50),
  textMatrix = textMatrix,
  setStdMargins = FALSE,
  cex.text = 0.5,
  cex.lab = 0.5,
  zlim = c(-1, 1),
  main = paste("Module-trait relationships")
)

geneNames <- data.frame(symbol = colnames(datExpr))
text <- unique(moduleColors)
for (i in 1:length(text)) {
  inModule <- is.finite(match(moduleColors, text[i]))
  modGenes <- geneNames[inModule, ]
  modGenes <- data.frame(modGenes)
  write.table(modGenes, paste(text[i], ".", "txt", sep = ""), sep = "\t", row.names = F, quote = F, col.names = F)
}

up_mod <- c("black", "brown")
down_mod <- c("salmon", "lightgreen")


####################################################
###################    xSum    #####################
####################################################
library(tidyverse)
library(data.table)
library(PharmacoGx)
library(parallel)
library(dplyr)
library(stringr)
library(tidyverse)
library(tibble)
library(clusterProfiler)
library(ggplot2)
library(ggrepel)
library(ggthemes)
library(gridExtra)

load("CMAP_gene_signatures.RData")
camp_sig <- CMAP.genePerturbations[, , c("tstat")] %>% data.frame()


camp_sig$ENTREZID <- do.call(rbind, strsplit(rownames(camp_sig), "\\."))[, 2]
SYMBOL <- bitr(camp_sig$ENTREZID, fromType = "ENTREZID", toType = "SYMBOL", OrgDb = "org.Hs.eg.db")
camp_sig <- merge(SYMBOL, camp_sig, by = "ENTREZID")
camp_sig <- column_to_rownames(camp_sig, var = "SYMBOL")
camp_sig <- camp_sig[, -1]




dd_up <- list()
for (i in list.files(pattern = "txt")) {
  if (is.element(gsub("\\.txt", "", i), up_mod)) {
    dd_up[[paste0("", gsub("\\.txt", "", i))]] <- fread(paste0("", i), data.table = F, header = F)
  } else {
    next
  }
}

dd_down <- list()
for (i in list.files(pattern = "txt")) {
  if (is.element(gsub("\\.txt", "", i), down_mod)) {
    dd_down[[paste0("", gsub("\\.txt", "", i))]] <- fread(paste0("", i), data.table = F, header = F)
  } else {
    next
  }
}

dd_up2 <- dd_up %>%
  unlist() %>%
  as.data.frame() %>%
  pull(.) %>%
  data.frame(id = .) %>%
  cbind(., fc = rep(1, length(dd_up)))
dd_down2 <- dd_down %>%
  unlist() %>%
  as.data.frame() %>%
  pull(.) %>%
  data.frame(id = .) %>%
  cbind(., fc = rep(-1, length(dd_down)))
dis_sig <- rbind(dd_up2, dd_down2)


camp_sig <- readRDS("camp_sig.rds")
source("xSum_Galaxy.R")

XLogFC <- eXtremeLogFC(camp_sig, N = 1000)
up_gene <- dis_sig$id[dis_sig$fc == 1]
dn_gene <- dis_sig$id[dis_sig$fc == -1]

xsum <- data.frame(score = XSum(XLogFC, up_gene, dn_gene))
xsum <- rownames_to_column(xsum, var = "id")
xsum <- xsum[order(xsum$score), ]

xsum$number <- 1:nrow(xsum)
write.csv(xsum, "xsum.csv")
select <- xsum[1:5, ]

ggplot(xsum, aes(number, score)) +
  geom_point(size = 3, color = "grey50") +
  geom_point(
    data = select, alpha = 1,
    size = 5, color = "#5ec7dd"
  ) +
  geom_label_repel(
    data = select, aes(label = id),
    color = "white",
    alpha = 1, point.padding = 1,
    size = 6.5, fill = "#009bc7",
    segment.size = 1, nudge_x = -0.5,
    segment.color = "grey50",
    direction = "x",
    hjust = 1
  ) +
  labs(y = "Xsum score", x = "Compounds") +
  theme_bw(base_rect_size = 1.5)












