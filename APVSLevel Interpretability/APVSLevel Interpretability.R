


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
p3 <- plotMsigNetwork_zg(gs_ovnet,
  markGroups = grps[1:8],
  nodeSF = 1.5,
  edgeSF = 0.001
)


source("plotMsigWordcloud_galaxy.R")
set.seed(21)
p4 <- plotMsigWordcloud_zg(msigGsc = siggs, groups = grps[1:8], type = "Name")
set.seed(21)
p5 <- plotMsigWordcloud_zg(msigGsc = siggs, groups = grps[1:8], type = "Short")


source("plotGeneStats_galaxy.R")
set.seed(21)
p6 <- plotGeneStats_zg(gStats,
  msigGsc = siggs,
  groups = grps[1:8],
  statName = "logFC"
)






####################################################
################    Immune    ######################
####################################################




ipsmap<- function (x) {
  if (x<=0) {
    ips<-0
  } else {
    if (x>=3) {
      ips<-10
    } else {
      ips<-round(x*10/3, digits=0)
    }
  }
  return(ips)
}



















