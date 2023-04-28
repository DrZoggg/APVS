


if(T){
  library(patchwork)
  library(ggalluvial)
  library(patchwork)
  library(reshape2)
  library(scales)
  library(ggpubr)
  library(pheatmap)
  library(hdf5r)
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
  library(risksetROC)
  library(survivalROC)
  library(survcomp)
  library(survminer)
  library(survival)
  library(ggkm)
  library(smoothHR)
  library(rms)
  library(ggrcs)
  library(qusage)
  library(spiralize)
  library(doParallel)
  library(pROC)
  library(vegan)
  library(tidyverse)
  library(dplyr)
  library(data.table)
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
##################   prognosis      ################
####################################################

pheno <- dd2[["Cohort1"]] %>%
  dplyr::rename("MACE" = time) %>%
  mutate(
    MACE_Censor = case_when(
      fusta %in% "TRUE" ~ 1, fusta %in% "FALSE" ~ 0
    ),
    MACE = MACE / 30
  )
BA12 <- dd1_CRDscore[["Cohort1"]]
rt <- data.frame(pheno[, c(5:6)], APVSLevel = BA12$APVSLevel) %>% 
  arrange(APVSLevel, desc(MACE))
res.cat <- surv_cutpoint(rt,
  time = "MACE", event = "MACE_Censor",
  variables = names(rt)[3:ncol(rt)]
) %>% surv_categorize()
cut_off <- surv_cutpoint(rt,
  time = "MACE", event = "MACE_Censor",
  variables = names(rt)[3:ncol(rt)]
)$cutpoint$cutpoint

group <- factor(res.cat$APVSLevel, levels = c("low", "high"))
my.surv <- Surv(rt$MACE, rt$MACE_Censor)
data.survdiff <- survdiff(my.surv ~ group)
p.val <- 1 - pchisq(data.survdiff$chisq, length(data.survdiff$n) - 1)
HR <- (data.survdiff$obs[2] / data.survdiff$exp[2]) / (data.survdiff$obs[1] / data.survdiff$exp[1])
up95 <- exp(log(HR) + qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
low95 <- exp(log(HR) - qnorm(0.975) * sqrt(1 / data.survdiff$exp[2] + 1 / data.survdiff$exp[1]))
HR <- paste("Hazard Ratio = ", round(HR, 2), sep = "")
CI <- paste("95% CI: ", paste(round(low95, 2), round(up95, 2), sep = " - "), sep = "")

KMplot <- ggplot(res.cat, aes(time = MACE, status = MACE_Censor, color = APVSLevel)) +
  geom_km(trans = "event") +
  geom_kmticks(trans = "event") +
  scale_color_manual(values = scales::alpha(c("#2C8EAE", "grey68"), 1)) +
  labs(y = "MACE", x = "Follow-up Time (month)") +
  ggtitle("APVSLevel") +
  annotate("text",
    x = 60, y = 0.04,
    label = paste(pval = ifelse(p.val < 0.0001, "Log-rank\np < 0.0001",
      paste("Log-rank\np = ", round(p.val, 5), sep = "")
    )),
    hjust = 0, size = 5.2, fontface = "bold.italic"
  ) +
coord_cartesian(clip = "on")




####################################################
################   fluctuation      ################
####################################################

GSE <- c("Cohort1", "Cohort2", "Cohort3", "Cohort4", "Cohort5", "Cohort6", "Cohort7")
dd_final <- list()
for (i in GSE) {
  BA12 <- dd1_APVSLevel[[i]]
  pheno <- dd2[[i]]
  data <- cbind(BA12, pheno)
  dd_final[[i]] <- data
}

ggdata <- data.frame()
ggdata_list <- list()
for (i in 1:length(dd_final)) {
  data <- dd_final[[i]] %>%
    group_by(group) %>%
    summarise(mean = mean(APVSLevel)) %>%
    as.data.frame() %>%
    cbind(., data.frame(batch = rep(names(dd_final)[i], nrow(.))))

  ggdata_list[[names(dd_final)[i]]] <- data
  ggdata <- rbind(ggdata, data)
}

ggdata_list <- lapply(ggdata_list, function(x) {
  data <- x
  data$group <- data$group %>% as.numeric()
  data <- data %>% arrange(group)
  return(data)
})

source('line_galaxy.R')
id = 'Cohort1'
p <- line_galaxy(data = ggdata_list[[id]],'#E69F00',2.5,5,id);p
ggsave(p,filename=paste0('ggline_',id,'.pdf'),width = 3.8,height= 2.2)
id = 'Cohort2'
p <- line_galaxy(data = ggdata_list[[id]],'#E69F00',2.5,5,id);p
ggsave(p,filename=paste0('ggline_',id,'.pdf'),width = 3.8,height= 2.2)
id = 'Cohort3'
p <- line_galaxy(data = ggdata_list[[id]],'#E69F00',2.5,5,id);p
ggsave(p,filename=paste0('ggline_',id,'.pdf'),width = 3.8,height= 2.2)
id = 'Cohort4'
p <- line_galaxy(data = ggdata_list[[id]],'#E69F00',2.5,5,id);p
ggsave(p,filename=paste0('ggline_',id,'.pdf'),width = 3.8,height= 2.2)
id = 'Cohort5'
p <- line_galaxy(data = ggdata_list[[id]],'#E69F00',2.5,5,id);p
ggsave(p,filename=paste0('ggline_',id,'.pdf'),width = 3.8,height= 2.2)
id = 'Cohort6'
p <- line_galaxy(data = ggdata_list[[id]],'#E69F00',2.5,5,id);p
ggsave(p,filename=paste0('ggline_',id,'.pdf'),width = 3.8,height= 2.2)
id = 'Cohort7'
p <- line_galaxy(data = ggdata_list[[id]],'#E69F00',2.5,5,id);p
ggsave(p,filename=paste0('ggline_',id,'.pdf'),width = 3.8,height= 2.2)



####################################################
##################   Severity      #################
####################################################


GSE <- "Cohort8"
dd_final <- list()
for (i in GSE) {
  BA12 <- dd1_APVSLevel[[i]]
  pheno <- dd2[[i]]
  data <- cbind(BA12, pheno)
  data$APVSLevel <- -data$APVSLevel
  data <- data %>% arrange(APVSLevel)
  data <- data %>% mutate(
    Degree = case_when(
      APVSLevel >= median(data$APVSLevel) ~ "Intensify",
      APVSLevel < median(data$APVSLevel) ~ "Alleviate"
    )
  )
  cut <- quantile(data$APVSLevel, 
  probs = seq(0,1, 0.1))
  data <- data %>% mutate(
    quantile = case_when(
      data$APVSLevel >= cut10[1] & data$APVSLevel < cut[2] ~ "0-10",
      data$APVSLevel >= cut10[2] & data$APVSLevel < cut[3] ~ "10-20",
      data$APVSLevel >= cut10[3] & data$APVSLevel < cut[4] ~ "20-30",
      data$APVSLevel >= cut10[4] & data$APVSLevel <= cut[5] ~ "30-40",
      data$APVSLevel >= cut10[5] & data$APVSLevel <= cut[6] ~ "40-50",
      data$APVSLevel >= cut10[6] & data$APVSLevel <= cut[7] ~ "50-60",
      data$APVSLevel >= cut10[7] & data$APVSLevel <= cut[8] ~ "60-70",
      data$APVSLevel >= cut10[8] & data$APVSLevel <= cut[9] ~ "70-80",
      data$APVSLevel >= cut10[9] & data$APVSLevel <= cut[10] ~ "80-90",
      data$APVSLevel >= cut10[10] & data$APVSLevel <= cut[11] ~ "90-100"
    )
  )
  dd_final[[i]] <- data
}

dt <- dd_final[[i]]
p <- ggviolin(dt,
  x = "quantile", y = "Disease index",
  fill = "quantile",
  alpha = .6,
  border = F,
  width = .95,
  add = c("median_mad"),
  add.params = list(color = "#250902", size = 0.4),
  draw_quantiles = .5,
  color = "white"
) +
  geom_point(
    size = 1, shape = 20,
    color = "grey10",
    position = position_jitterdodge(),
    aes(fill = "quantile"), alpha = 1
  ) +
  scale_fill_manual(values = alpha(colorRampPalette(paletteer_d("ggsci::blue_grey_material"),
    alpha = TRUE
  )(n = 11), 1)) +
  scale_color_manual(values = alpha(colorRampPalette(paletteer_c("grDevices::Set 3", 11),
    alpha = TRUE
  )(n = 11), 1)) +
  xlab("APVSLevel") +
  theme_classic() +
  scale_y_continuous(limits = c(-0, 120), expand = c(0.02, 0)) +
stat_compare_means(size = 4)




####################################################
#############    Extrapolability    ################
####################################################


GSE <- c("Cohort1","Cohort2","Cohort3","Cohort4","Cohort5","Cohort6","Cohort7","Cohort8","Cohort9","Cohort10")

dd_final <- list()
for (i in GSE) {
  BA12 <- dd1_APVSLevel[[i]]
  pheno <- dd2[[i]]
  data <- cbind(BA12, pheno)
  data <- data %>% arrange(APVSLevel)
  data$batch <- i
  data <- data %>% dplyr::select(APVSLevel, batch)
  data <- data %>% rownames_to_column("id")
  dd_final[[i]] <- data
}

mpg <- data.frame()
for (i in names(dd_final)) {
  rt <- dd_final[[i]]
  mpg <- rbind(mpg, rt)
}

cut <- quantile(mpg$APVSLevel, 
       probs = seq(0,1,0.1))
mpg <- mpg %>% mutate(
  Degree = case_when(
    APVSLevel >= median(mpg$APVSLevel) ~ "APVSLevel-High",
    APVSLevel < median(mpg$APVSLevel) ~ "APVSLevel-Low"
  ),
  quantile = case_when(
    mpg$APVSLevel >= cut10[1] & mpg$APVSLevel < cut[2] ~ "0-10",
    mpg$APVSLevel >= cut10[2] & mpg$APVSLevel < cut[3] ~ "10-20",
    mpg$APVSLevel >= cut10[3] & mpg$APVSLevel < cut[4] ~ "20-30",
    mpg$APVSLevel >= cut10[4] & mpg$APVSLevel <= cut[5] ~ "30-40",
    mpg$APVSLevel >= cut10[5] & mpg$APVSLevel <= cut[6] ~ "40-50",
    mpg$APVSLevel >= cut10[6] & mpg$APVSLevel <= cut[7] ~ "50-60",
    mpg$APVSLevel >= cut10[7] & mpg$APVSLevel <= cut[8] ~ "60-70",
    mpg$APVSLevel >= cut10[8] & mpg$APVSLevel <= cut[9] ~ "70-80",
    mpg$APVSLevel >= cut10[9] & mpg$APVSLevel <= cut[10] ~ "80-90",
    mpg$APVSLevel >= cut10[10] & mpg$APVSLevel <= cut[11] ~ "90-100"
  )
)

p <- ggplot(mpg, aes(batch)) +
  geom_bar(aes(fill = Degree),
    position = "fill",
    width = .9, color = "black", size = .3
  ) +
  coord_polar(start = 0) +
  labs(x = "batch", y = "Percentage");P

pl <- list()
for (i in names(dd_final)) {
  dt <- dd_final[[i]] %>%
    arrange(APVSLevel) %>%
    mutate(num = 1:nrow(.))
  pl[[i]] <- ggplot(dt, aes(x = num, y = APVSLevel)) +
    geom_point(
      size = 1.7,
      shape = 21, alpha = .8,
      fill = mycolor[[i]],
      color = mycolor[[i]]
    ) +
    geom_smooth(
      method = "gam",
      size = 1.3,
      alpha = 1,
      color = mycolor[[i]],
      se = F
    ) +
    labs(title = i) +
    theme_cleantable() +
    annotate("segment",
      x = quantile(dt$num)[2], xend = quantile(dt$num)[4],
      y = mean(dt$APVSLevel), yend = mean(dt$APVSLevel),
      size = 1, alpha = 1, color = "black"
    ) +
    annotate("segment",
      x = 1, xend = max(dt$num), lty = 2,
      y = median(mpg$APVSLevel), yend = median(mpg$APVSLevel),
      size = 1.5, alpha = 1, color = "grey70"
    )
}

p_combined <- ggarrange(
  plotlist = pl,
  nrow = 1, ncol = length(pl),
  align = "v",
  label.y = 0,
  heights = c(2, 2),
  common.legend = F
)



####################################################
########    Biological implications    #############
####################################################
library(msigdb)
library(msigdbr)
library(GSEABase)
library(vissE)
library(igraph)
library(qusage)


GSE <- c("Cohort1", "Cohort2", "Cohort3", "Cohort4", "Cohort5", "Cohort6", "Cohort7", "Cohort8", "Cohort9", "Cohort10")
dd_final2 <- list()
for (i in GSE) {
  BA12 <- dd1[[i]] %>%
    t() %>%
    as.data.frame()
  BA12 <- cbind(BA12, dd1_APVSLevel[[i]])
  data <- BA12
  data <- data %>% arrange(APVSLevel)
  data$batch <- i
  dd_final2[[i]] <- data
}

dd_genelist <- lapply(dd_final2, function(x) {
  x <- colnames(x)
  return(x)
})
dd_final3 <- lapply(dd_final2, function(x) {
  data <- x[, Reduce(intersect, dd_genelist)]
})

mpg <- data.frame()
for (i in names(dd_final3)) {
  rt <- dd_final3[[i]]
  mpg <- rbind(mpg, rt)
}

mpg <- mpg %>% mutate(
  Degree = case_when(
    APVSLevel >= median(mpg$APVSLevel) ~ "Intensify",
    APVSLevel < median(mpg$APVSLevel) ~ "Alleviate"
  )
)

phenotypes <- mpg %>%
  dplyr::select(c(Degree)) %>%
  pull(Degree)
phenotypes <- factor(phenotypes, levels = c("Intensify", "Alleviate"))
phenotypes
exprsdata <- mpg %>%
  dplyr::select(-c("APVSLevel", "batch", "Degree")) %>%
  t() %>%
  as.matrix()

n.sub <- length(table(phenotypes))
n.sub.label <- unique(phenotypes)
treat_list <- ctrl_list <- degs.list <- list()

for (i in 1:n.sub) {
  cat(paste0(n.sub.label[i], " vs. Others go!\n"))
  treat_list[[i]] <- phenotypes[which(phenotypes == n.sub.label[i])]
  ctrl_list[[i]] <- phenotypes[-which(phenotypes == n.sub.label[i])]

  treat <- treat_list[[i]]
  ctrl <- ctrl_list[[i]]

  mydata <- exprsdata[, c(treat, ctrl)]
  group_list <- factor(c(
    rep("treat", length(treat)),
    rep("ctrl", length(ctrl))
  ))
  design <- model.matrix(~ 0 + group_list)
  colnames(design) <- levels(group_list)
  rownames(design) <- colnames(mydata)
  contrst <- makeContrasts("treat-ctrl", levels = design)
  fit <- lmFit(mydata, design = design)
  fit2 <- contrasts.fit(fit, contrasts = contrst)
  fit2 <- eBayes(fit2)
  alldiff <- topTable(fit2, coef = 1, Inf)

  degs.list[[n.sub.label[i]]] <- as.data.frame(na.omit(alldiff))
  cat("\n")
}

sig.markers <- degs.list[[1]] %>% filter(P.Value < 0.01 & abs(logFC) > 1)

load("msigdb_hs.rda")
msigdb_hs
genesigs <- c(
  subsetCollection(msigdb_hs)
)
genesigs <- GeneSetCollection(genesigs)
hallmark_geneset_ls <- geneIds(genesigs)


df.data <- exprsdata
labels <- phenotypes
contrast <- paste("Intensify", "Alleviate", sep = "-")
set.seed(123)
res_qusage_hallmark <- qusage(
  eset = df.data,
  labels = labels,
  contrast = contrast,
  geneSets = hallmark_geneset_ls
)


res_table <- qsTable(res_qusage_hallmark, number = 50000, sort.by = "logFC")
res_table <- res_table %>% mutate(
  Direction = case_when(
    res_table$log.fold.change >= 0 ~ "Up",
    res_table$log.fold.change < 0 ~ "Down",
  )
)
res_table_sig <- res_table %>% filter(FDR < 0.05 & abs(log.fold.change) > 0.5)
siggs <- genesigs[res_table_sig$pathway.name]


set.seed(21)
gs_ovlap <- computeMsigOverlap(siggs, thresh = 0.1, measure = "jaccard")
gs_ovnet <- computeMsigNetwork(gs_ovlap, msigGsc = siggs)
set.seed(211)
p1 <- plotMsigNetwork(gs_ovnet, nodeSF = 1.5, edgeSF = 0.1)


gsStats <- res_table_sig$log.fold.change
names(gsStats) <- res_table_sig$pathway.name
gStats <- degs.list[[1]]$logFC
names(gStats) <- degs.list[[1]] %>% row.names()
set.seed(211)
p2 <- plotMsigNetwork(gs_ovnet, genesetStat = gsStats, nodeSF = 1.5, edgeSF = 0.1)


grps <- cluster_walktrap(gs_ovnet)
grps <- groups(grps)
grps <- grps[sapply(grps, length) > 10]
grp_size <- sapply(grps, length)
grp_stat <- sapply(grps, function(x) median(abs(gsStats[x])))
grp_pr <- rank(grp_size) * rank(grp_stat)
grps <- grps[order(grp_pr, decreasing = TRUE)]
source("plotMsigNetwork_galaxy.R")
p3 <- plotMsigNetwork_galaxy(gs_ovnet,
  markGroups = grps[1:8],
  nodeSF = 1.5,
  edgeSF = 0.001
)


source("plotMsigWordcloud_galaxy.R")
set.seed(21)
p4 <- plotMsigWordcloud_galaxy(msigGsc = siggs, groups = grps[1:8], type = "Name")
set.seed(21)
p5 <- plotMsigWordcloud_galaxy(msigGsc = siggs, groups = grps[1:8], type = "Short")


source("plotGeneStats_galaxy.R")
set.seed(21)
p6 <- plotGeneStats_galaxy(gStats,
  msigGsc = siggs,
  groups = grps[1:8],
  statName = "logFC"
)


source('plotMsigPPI_galaxy.R')
getIMEXVersions()
ppi = getIMEX(org = 'hs', inferred = TRUE)
set.seed(21)
p7 <- plotMsigPPI_galaxy(
  ppidf = ppi,
  msigGsc = siggs,
  groups = grps[1:8],
  geneStat = gStats,
  threshFrequency = 0.0, 
  threshConfidence = 0.0, 
  threshStatistic = 0.0, 
  threshUseAbsolute = TRUE 
)



####################################################
################    Immune    ######################
####################################################
GSE <- c("Cohort1", "Cohort2", "Cohort3", "Cohort4", "Cohort5", "Cohort6", "Cohort7", "Cohort8", "Cohort9", "Cohort10")
dd_final2 <- list()
for (i in GSE) {
  BA12 <- dd1[[i]] %>%
    t() %>%
    as.data.frame()
  BA12 <- cbind(BA12, dd1_APVSLevel[[i]])
  data <- BA12
  data <- data %>% arrange(APVSLevel)
  data$batch <- i
  dd_final2[[i]] <- data
}

dd_genelist <- lapply(dd_final2, function(x) {
  x <- colnames(x)
  return(x)
})
dd_final3 <- lapply(dd_final2, function(x) {
  data <- x[, Reduce(intersect, dd_genelist)]
})

mpg <- data.frame()
for (i in names(dd_final3)) {
  rt <- dd_final3[[i]]
  mpg <- rbind(mpg, rt)
}

mpg <- mpg %>% mutate(
  Degree = case_when(
    APVSLevel >= median(mpg$APVSLevel) ~ "Intensify",
    APVSLevel < median(mpg$APVSLevel) ~ "Alleviate"
  )
)

phenotypes <- mpg %>%
  dplyr::select(c(Degree)) %>%
  pull(Degree)
phenotypes <- factor(phenotypes, levels = c("Intensify", "Alleviate"))
phenotypes
exprsdata <- mpg %>%
  dplyr::select(-c("APVSLevel", "batch", "Degree")) %>%
  t() %>%
  as.matrix()

tumsam <- colnames(exprsdata)
gene_expression <- exprsdata # 对TPM值做对数转化
sample_names <- tumsam



load("IPS ref.rda")
head(IPSG)
unique_ips_genes <- as.vector(unique(IPSG$NAME))
IPS <- MHC <- CP <- EC <- SC <- AZ <- NULL
GVEC <- row.names(gene_expression)
VEC <- as.vector(IPSG$GENE)

ind <- which(is.na(match(VEC, GVEC)))
MISSING_GENES <- VEC[ind]
MISSING_GENES

which(GVEC %in% VEC)
dat <- IPSG[ind, ]
if (length(MISSING_GENES) > 0) {
  message(paste0("--differently named or missing genes: ", paste(MISSING_GENES, collapse = ", ")))
  print(IPSG[ind, ])
  message("please check data and make sure all genes matches!")
} else {
  message("--all genes matched!")
}
# --differently named or missing genes: IARS, CD300E, DARS, SIGLEC14, FCRL6, SIK1, FCGR3A, CCL3L1
# GENE    NAME CLASS WEIGHT
# 34      IARS Act CD4    EC      1
# 74    CD300E Tem CD4    EC      1
# 75      DARS Tem CD4    EC      1
# 88  SIGLEC14 Tem CD4    EC      1
# 109    FCRL6 Tem CD8    EC      1
# 121     SIK1 Tem CD8    EC      1
# 130   FCGR3A    MDSC    SC     -1
# 143   CCL3L1    Treg    SC     -1
# please check data and make sure all genes matches!


unique_ips_genes <- as.vector(unique(IPSG$NAME))
IPS <- MHC <- CP <- EC <- SC <- AZ <- NULL
GVEC <- row.names(gene_expression)
VEC <- as.vector(IPSG$GENE)
ind <- which(is.na(match(VEC, GVEC)))
MISSING_GENES <- VEC[ind]
MISSING_GENES



my_palette <- colorRampPalette(c("#0F7BA2", "white", "#DD5129"))(n = 1000)
mapcolors <- function(x) {
  za <- NULL
  if (x >= 3) {
    za <- 1000
  } else {
    if (x <= -3) {
      za <- 1
    } else {
      za <- round(166.5 * x + 500.5, digits = 0)
    }
  }
  return(my_palette[za])
}


my_palette2 <- colorRampPalette(c("white", "#90e0ef"))(n = 1000)
mapbw <- function(x) {
  za2 <- NULL
  if (x >= 2) {
    za2 <- 1000
  } else {
    if (x <= -2) {
      za2 <- 1
    } else {
      za2 <- round(249.75 * x + 500.5, digits = 0)
    }
  }
  return(my_palette2[za2])
}


ipsmap <- function(x) {
  if (x <= 0) {
    ips <- 0
  } else {
    if (x >= 3) {
      ips <- 10
    } else {
      ips <- round(x * 10 / 3, digits = 0)
    }
  }
  return(ips)
}





library(gridExtra)
outTab <- NULL
for (i in 1:length(sample_names)) { 
  GE <- gene_expression[[i]]
  mGE <- mean(GE)
  sGE <- sd(GE)
  Z1 <- (gene_expression[as.vector(IPSG$GENE),i]-mGE)/sGE
  W1 <- IPSG$WEIGHT
  WEIGHT <- MIG <- NULL
  k <- 1
  for (gen in unique_ips_genes) {
    MIG[k] <- mean(Z1[which (as.vector(IPSG$NAME)==gen)],na.rm=TRUE)
    WEIGHT[k] <- mean(W1[which (as.vector(IPSG$NAME)==gen)])
    k <- k + 1
  }
  WG <- MIG * WEIGHT
  MHC[i] <- mean(WG[1:10])
  CP[i] <- mean(WG[11:20])
  EC[i] <- mean(WG[21:24])
  SC[i] <- mean(WG[25:26])
  AZ[i] <- MHC[i]+CP[i]+EC[i]-SC[i]
  IPS[i] <- ipsmap(AZ[i])
  
  tmp <- as.data.frame(t(data.frame(WG))); colnames(tmp) <- unique_ips_genes; rownames(tmp) <- sample_names[i]
  outTab <- rbind.data.frame(outTab,tmp)
  
  data_a <- data.frame (start = c(0,2.5,5,7.5,10,15,seq(20,39),0,10,20,30), end = c(2.5,5,7.5,10,15,seq(20,40),10,20,30,40), y1=c(rep(2.6,26),rep(0.4,4)),y2=c(rep(5.6,26),rep(2.2,4)),z=c(MIG[c(21:26,11:20,1:10)],EC[i],SC[i],CP[i],MHC[i]),vcol=c(unlist(lapply(MIG[c(21:26,11:20,1:10)],mapcolors)), unlist(lapply(c(EC[i],SC[i],CP[i],MHC[i]),mapbw))), label = c(unique_ips_genes[c(21:26,11:20,1:10)],"EC","SC","CP","MHC"))
  data_a$label <- factor(data_a$label, levels=unique(data_a$label))
  plot_a1 <- ggplot() + 
    geom_rect(data=data_a, mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +
    coord_polar()  + 
    scale_y_continuous(limits = c(0, 6)) +
    scale_fill_manual(values =as.vector(data_a$vcol),guide=FALSE) +
    theme_bw() + 
    theme(panel.spacing = unit(0, 'mm'), 
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          axis.line = element_line(colour = "white"),
          axis.text=element_blank(), 
          axis.ticks= element_blank()) + 
    geom_text(aes(x=5, y=1.3, label="EC"), size=4) +
    geom_text(aes(x=15, y=1.3, label="SC"), size=4) + 
    geom_text(aes(x=25, y=1.3, label="CP"), size=4) + 
    geom_text(aes(x=35, y=1.3, label="MHC"), size=4)
  plot_a2 <- plot_a1 +
    geom_text(aes(x=1.25, y=4.1, label="+ Act CD4"), angle=78.75, size=4) +
    geom_text(aes(x=3.75, y=4.1, label="+ Act CD8"),angle=56.25, size=4) +
    geom_text(aes(x=6.25, y=4.1, label="+ Tem CD4"), angle=33.75,size=4) +
    geom_text(aes(x=8.75, y=4.1, label="+ Tem CD8"), angle=11.25,size=4) +
    geom_text(aes(x=17.5, y=4.1, label="- MDSC"), angle=-67.5,size=4) +
    geom_text(aes(x=12.5, y=4.1, label="- Treg"), angle=-22.5,size=4)
  plot_a3 <- plot_a2 +
    geom_text(aes(x=20.5, y=4.1, label="PD-1 -"), angle=85.5, size=4) +
    geom_text(aes(x=21.5, y=4.1, label="CTLA4 -"), angle=76.5, size=4) +
    geom_text(aes(x=22.5, y=4.1, label="LAG3 -"), angle=67.5, size=4) +
    geom_text(aes(x=23.5, y=4.1, label="TIGIT -"), angle=58.5, size=4) +
    geom_text(aes(x=24.5, y=4.1, label="TIM3 -"), angle=49.5, size=4) +
    geom_text(aes(x=25.5, y=4.1, label="PD-L1 -"), angle=40.5, size=4) +
    geom_text(aes(x=26.5, y=4.1, label="PD-L2 -"), angle=31.5, size=4) +
    geom_text(aes(x=27.5, y=4.1, label="CD27 +"), angle=22.5, size=4) +
    geom_text(aes(x=28.5, y=4.1, label="ICOS +"), angle=13.5, size=4) +
    geom_text(aes(x=29.5, y=4.1, label="IDO1 -"), angle=4.5, size=4)
  plot_a4 <- plot_a3 +
    geom_text(aes(x=30.5, y=4.1, label="B2M +"), angle=-4.5, size=4) +
    geom_text(aes(x=31.5, y=4.1, label="TAP1 +"), angle=-13.5, size=4) +
    geom_text(aes(x=32.5, y=4.1, label="TAP2 +"), angle=-22.5, size=4) +
    geom_text(aes(x=33.5, y=4.1, label="HLA-A +"), angle=-31.5, size=4) +
    geom_text(aes(x=34.5, y=4.1, label="HLA-B +"), angle=-40.5, size=4) +
    geom_text(aes(x=35.5, y=4.1, label="HLA-C +"), angle=-49.5, size=4) +
    geom_text(aes(x=36.5, y=4.1, label="HLA-DPA1 +"), angle=-58.5, size=4) +
    geom_text(aes(x=37.5, y=4.1, label="HLA-DPB1 +"), angle=-67.5, size=4) +
    geom_text(aes(x=38.5, y=4.1, label="HLA-E +"), angle=-76.5, size=4) +
    geom_text(aes(x=39.5, y=4.1, label="HLA-F +"), angle=-85.5, size=4)
  plot_a5 <- plot_a4 +
    geom_text(aes(x=0, y=6, label=paste("Immunophenoscore: ",IPS[i],sep="")), angle=0,size=6,vjust=-0.5) + 
    theme(axis.title=element_blank())
  plot_a <- plot_a5 + theme(plot.margin=unit(c(0,0,0,0),"mm")) +
    geom_text(vjust=1.15,hjust=0,aes(x=25.5, y=6,label="\n\n\n\n   MHC: Antigen Processing                                 EC: Effector Cells\n   CP: Checkpoints | Immunomodulators              SC: Suppressor Cells\n\n", hjust = 0), size=4)
  
  data_b <- data.frame (start = rep(0,23), 
                        end = rep(0.7,23), 
                        y1=seq(0,22,by=1), 
                        y2=seq(1,23,by=1),
                        z=seq(-3,3,by=6/22),
                        vcol=c(unlist(lapply(seq(-3,3,by=6/22),mapcolors))), 
                        label = LETTERS[1:23])
  data_b_ticks <- data.frame(x = rep(1.2, 7), 
                             value = seq(-3,3, by=1), 
                             y = seq(0,6, by=1)*(22/6) +0.5)
  legendtheme <- theme(plot.margin = unit(c(2,0,2,0),"inch"), 
                       panel.spacing = unit(0,"null"), 
                       panel.grid.major = element_blank(),
                       panel.grid.minor = element_blank(),
                       panel.border = element_blank(),
                       panel.background = element_blank(), 
                       axis.line = element_line(colour = "white"), 
                       axis.text=element_blank(), 
                       axis.ticks= element_blank(), 
                       axis.title.x=element_blank())
  plot_b <- ggplot(hjust=0) + 
    geom_rect(data=data_b, 
              mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label), size=0.5,color="black", alpha=1) +
    scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) +
    scale_fill_manual(values =as.vector(data_b$vcol),guide=FALSE) +
    geom_text(data=data_b_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +
    theme_bw() +
    legendtheme + 
    ylab("Sample-wise (averaged) z-score")
  
  data_c <- data.frame(start = rep(0,23),
                       end = rep(0.7,23),
                       y1=seq(0,22,by=1),
                       y2=seq(1,23,by=1),
                       z=seq(-2,2,by=4/22),
                       vcol=c(unlist(lapply(seq(-2,2,by=4/22),mapbw))), 
                       label = LETTERS[1:23])
  data_c_ticks <- data.frame(x = rep(1.2, 5), 
                             value = seq(-2,2, by=1),
                             y = seq(0,4, by=1)*(22/4) +0.5)
  plot_c <- ggplot() + 
    geom_rect(data=data_c, 
              mapping=aes(xmin=start, xmax=end, ymin=y1, ymax=y2, fill=label),size=0.5,color="black", alpha=1) + 
    scale_x_continuous(limits = c(0, 1.5),expand = c(0,0)) + 
    scale_fill_manual(values =as.vector(data_c$vcol),guide=FALSE) + 
    geom_text(data=data_c_ticks, aes(x=x, y=y, label=value),hjust="inward", size=4) +
    theme_bw() + 
    legendtheme + 
    ylab("Weighted z-score")
  
  file_name<-paste("IPS_",sample_names[i],".pdf",sep="")
  pdf(file_name, width=6, height=5.5)
  grid.arrange(plot_a,plot_b,plot_c, ncol=3, widths=c(0.8,0.1,0.1))
  invisible(dev.off())
  
  message(paste0("--analysis of ", tumsam[i]," done..."))
}


DF <- data.frame(
  SAMPLE = sample_names,
  MHC = MHC,
  EC = EC,
  SC = SC,
  CP = CP,
  AZ = AZ,
  IPS = IPS,
  stringsAsFactors = F
)
write.table(outTab, file = "output_zscore.txt", row.names = TRUE, col.names = NA, quote = FALSE, sep = "\t")
write.table(DF, file = "output_IPS.txt", row.names = FALSE, col.names = TRUE, quote = FALSE, sep = "\t")


