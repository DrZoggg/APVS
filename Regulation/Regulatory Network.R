

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
}



mod_data <- fread("expr.csv", header = T) %>%
  filter(V1 %in% DCPGs) %>%
  column_to_rownames("V1") %>%
  t() %>%
  merge(., read.csv("clin.csv", row.names = 1, header = T) %>% select(group), by.x = 0, by.y = 0) %>%
  column_to_rownames("Row.names")


####################################################
##################   Hallmark      ################
####################################################

expr = mod_data%>%select(-c('group'))
load("hallmark.gs.RData")
gsva_es <- gsva(as.matrix(expr), gs)
row.names(gsva_es) <- str_replace(row.names(gsva_es), "HALLMARK_", "")
write.csv(gsva_es, "gsva_output_hall.csv", quote = F)


gsva_es <- read.csv("gsva_output_hall.csv",row.names = 1,check.names = F)
group_list <- mod_data %>% select(group) 
group_list$group <- factor(group_list$group,levels = c("STEMI","CCS"))
gsva_es = gsva_es[,rownames(group_list)]


design <- model.matrix(~ 0 + group_list$group)
colnames(design) <- levels(group_list$group)
rownames(design) <- colnames(gsva_es)
fit <- lmFit(gsva_es, design)
fit2 <- eBayes(fit)
alldiff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
write.csv(alldiff, "gsva_limma_hallmark.csv", quote = F)


gsva_es <- read.csv("gsva_output_hall.csv",
                    row.names = 1,check.names = F)
sigdiff <- subset(alldiff, abs(t) >= 5 &  
                    adj.P.Val < 0.05)
sig_gene <- rownames(sigdiff) 
write.table(sig_gene, file = 'sig_HALLMARK.txt')   
exp_key = gsva_es[sig_gene ,] 
write.csv(exp_key,'gsva_HALLMARK_key.csv')



#########==========heatmap
hmExp = read.csv( "gsva_HALLMARK_key.csv",row.names = 1)
plotdata <- as.matrix(hmExp)
targets = read.csv('targets.csv') 
targets$order = targets$group
targets$order[which(targets$order=='Stable CCS')] = 1
targets$order[which(targets$order=='1st day of SETMI')] = 2
targets$order[which(targets$order=='4-6 days of MI')] = 3
targets$order[which(targets$order=='1 month after MI')] = 4
targets$order[which(targets$order=='6 months after MI')] = 5
cc <- targets
cc <- cc %>%arrange(order)
cc$group = factor(cc$group, levels = c( 'Stable CCS','1st day of SETMI','4-6 days of SETMI','1 month after SETMI','6 months after SETMI'
))
plotdata = plotdata[, cc$id]

source('standarize.fun.R')
plotdata_mean = cbind(t(standarize.fun(plotdata)),data.frame(group = cc$group)) %>%
  group_by(group) %>% summarise_all(mean) %>% column_to_rownames('group')  %>% t() %>% as.data.frame()
Top_mean = HeatmapAnnotation(RISK = unique(cc$group),annotation_legend_param=list(labels_gp = gpar(fontsize = 10),
                                                          title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                          ncol=1),border = T,gap = unit(20.5, "cm"),
                             col=list(RISK = c( 'Stable CAD' = '#32A087',
                                                '1st day of MI' = '#DD492E',
                                                '4-6 days of MI' = "#40548A",
                                                '1 month after MI' = '#774f38',
                                                '6 months after MI' = '#EC7D21')),
                             show_annotation_name = F,annotation_name_side="left",
                             annotation_name_gp = gpar(fontsize = 10),simple_anno_size_adjust = T,show_legend = T,annotation_width = unit(0.5, "cm")
)
Left = rowAnnotation( Level = rep('HALLMARK',length( row.names(plotdata)))  ,
                      annotation_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"),ncol=1),
                      border = T,col=list(Level =  c( 'HALLMARK'=pal_npg("nrc", alpha = 0.6)(9)[3])),
                      show_annotation_name = F,annotation_name_gp = gpar(fontsize = 10),
                      simple_anno_size_adjust = T,show_legend = T,
                      annotation_width = unit(0.1, "cm")
)
hm1_mean = Heatmap(as.matrix(plotdata_mean),
                   name=' HALLMARK ',
                   left_annotation = Left,
                   cluster_rows = T,
                   cluster_columns = F,
                   color_space = "RGB",
                   col = colorRamp2(seq(-1.2, 1.2, length = 3), c("#00beef", 'white','#ea7d93')),
                   border = T,
                   border_gp = gpar(col = "black",lty = 1,lwd=1.5, linejoin = 'round' ),
                   rect_gp = gpar(col = "white",lty = 1,lwd=2, linejoin = 'round' ),
                   row_dend_side = c("left"),
                   row_dend_width = unit(1, "mm"),
                   row_order=NULL,
                   column_order=NULL,
                   column_names_rot = 45,
                   row_names_side = 'right',
                   show_column_names = F,
                   show_row_names = T,  
                   row_names_gp = gpar(col = 'grey30',fontsize = 7),
                   row_names_max_width = unit(12, "cm"),
                   column_labels = colnames(plotdata),
                   column_split =  unique(cc$group) ,
                   row_gap = unit(c(2), "mm"),
                   column_gap = unit(c(2), "mm"),
                   column_title = NULL,
                   column_title_gp = gpar(fill = "red", col = "white", border = "blue",fontface = "bold",fontsize = 9),
                   row_title_gp = gpar(fontsize=12),
                   show_heatmap_legend = T,
                   heatmap_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"))) 





####################################################
#####################   APVS      ##################
####################################################


#########==========heatmap
exp = mod_data%>%select(-c('group'))
keygene = read.table("APVS gene.txt") %>% .[, 1]
exp_key = exp[keygene$V1, ]
plotdata <- as.matrix(exp_key)
plotdata = plotdata[, cc$id]

plotdata_mean = cbind(t(standarize.fun(plotdata)),data.frame(group = cc$group)) %>%
  group_by(group) %>% 
  summarise_all(mean) %>% column_to_rownames('group')  %>% t() %>% as.data.frame()

Top_mean = HeatmapAnnotation(RISK = unique(cc$group),
                             annotation_legend_param=list(labels_gp = gpar(fontsize = 10),
                                                          title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                          ncol=1),
                             border = T,
                             gap = unit(30.5, "cm"),
                             col = list(RISK = my_palette),
                             show_annotation_name = F,
                             annotation_name_side="left",
                             annotation_name_gp = gpar(fontsize = 10),
                             simple_anno_size_adjust = T, 
                             show_legend = T,
                             annotation_width = unit(0.5, "cm")
)
Left = rowAnnotation( Level = rep('APVS', length( row.names(plotdata)))  ,
                      annotation_legend_param=list(labels_gp = gpar(fontsize = 10),
                                                   title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                   ncol=1),border = T, col=list(Level =   c( 'APVS' =  pal_npg("nrc", alpha = 0.6)(9)[1]   )           ),
                      show_annotation_name = F,annotation_name_gp = gpar(fontsize = 10),simple_anno_size_adjust = T,  
                      show_legend = T,annotation_width = unit(0.2, "cm")
)

hm2_mean = Heatmap(as.matrix(plotdata_mean),
                   name='APVS',
                   top_annotation = Top_mean,
                   left_annotation = Left,
                   cluster_rows = T,
                   cluster_columns = F,
                   color_space = "RGB",
                   col = colorRamp2(seq(-1.5, 1.5, length = 3), c("#84D2EE",'white','#FF7867')),
                   border = T,
                   border_gp = gpar(col = "black",lty = 1,lwd=1.5, linejoin = 'round' ),
                   rect_gp = gpar(col = "white",lty = 1,lwd=2, linejoin = 'round' ),
                   row_dend_side = c("left"),
                   row_dend_width = unit(1, "mm"),
                   row_order=NULL,
                   column_order=NULL,
                   column_names_rot = 45,
                   row_names_side = 'right',
                   show_column_names = F,
                   show_row_names = T,  
                   row_names_gp = gpar(col = 'grey30',fontsize = 7),
                   row_names_max_width = unit(12, "cm"),
                   column_labels = colnames(plotdata),
                   column_split =  unique(cc$group) ,
                   row_gap = unit(c(2), "mm"),
                   column_gap = unit(c(2), "mm"),
                   column_title = NULL,
                   column_title_gp = gpar(fill = "red", col = "white", border = "blue",fontface = "bold",fontsize = 9),
                   row_title_gp = gpar(fontsize=12),
                   show_heatmap_legend = T,
                   heatmap_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"))) 



####################################################
#####################   Cell      ##################
####################################################

expr = mod_data%>%select(-c('group'))
xCell = xCellAnalysis(t(scale(t(exp))),rnaseq=FALSE,scale =F,parallel.sz = 8) 
write.csv(xCell,file = "xCell.csv",row.names = T,quote = F)

cell_type = data.frame(
  cell_type = c('HSC','CLP','CMP','GMP','Megakaryocytes','Erythrocytes','Platelets',
                'pro B-cells','naive B-cells','B-cells','Class-switched memory B-cells','Memory B-cells','Plasma cells','Th1 cells','Th2 cells','Tregs' , 'CD4+ naive T-cells', 'CD4+ memory T-cells','CD4+ Tcm', 'CD4+ Tem', 'CD8+ naive T-cells','CD8+ T-cells','CD8+ Tcm','CD8+ Tem','NK cells','NKT' ,
                'Eosinophils','Neutrophils','Basophils','pDC','cDC','iDC','aDC', 'DC','Mast cells','Monocytes', 'Macrophages','Macrophages M1','Macrophages M2'),     
  group = c(rep('HSC',7),rep('Lymphoid',19),rep('Myeloid',13)    )
)
cell_type$cell_type %in% rownames(xCell)

gsva_es <- read.csv("xCell.csV",header = T,row.names = 1)
gsva_es = gsva_es[,rownames(group_list)]

design <- model.matrix(~ 0 + group_list$group)
colnames(design) <- levels(group_list$group)
rownames(design) <- colnames(gsva_es)
contrst <- makeContrasts('STEMI-CCS',levels = design)


fit <- lmFit(gsva_es, design)
fit2 <- contrasts.fit(fit,contrasts = contrst)
fit2 <- eBayes(fit2)
alldiff <- topTable(fit2,coef = 1,Inf)
head(alldiff)
sigdiff <- subset(alldiff,
                  abs(t) > 5 & 
                    adj.P.Val < 0.05)
sig_gene <- rownames(sigdiff)
write.table(sig_gene,'sig_xcell.txt' )



#########==========heatmap
exp = read.csv("xCell.csv",header = T,row.names = 1)
sig_gene = read.table('sig_xcell.txt')[,1]
exp_key = exp[sig_gene ,]
plotdata <- as.matrix(exp_key)

plotdata_mean = cbind(t(standarize.fun(plotdata)),data.frame(group = cc$group)) %>%group_by(group) %>% 
  summarise_all(mean) %>% column_to_rownames('group') %>% t() %>% as.data.frame()

Left = rowAnnotation( Level = rep('Cell', length( row.names(plotdata)))  ,
                      annotation_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"),
                      ncol=1),border = T, col=list(Level =   c( 'Cell' =  pal_npg("nrc", alpha = 0.6)(9)[4])),
                      show_annotation_name = F,annotation_name_gp = gpar(fontsize = 10),
                      simple_anno_size_adjust = T,show_legend = T,annotation_width = unit(0.2, "cm")
)


hm3_mean = Heatmap(as.matrix(plotdata_mean),
                   name='Cell',
                   left_annotation = Left,
                   cluster_rows = T,
                   cluster_columns = F,
                   color_space = "RGB",
                   col = colorRamp2(seq(-1.5, 1.5, length = 3), c("#88FFF7",'white','#7C83FD')),
                   border = T,
                   border_gp = gpar(col = "black",lty = 1,lwd=1.5, linejoin = 'round' ),
                   rect_gp = gpar(col = "white",lty = 1,lwd=2, linejoin = 'round' ),
                   row_dend_side = c("left"),
                   row_dend_width = unit(1, "mm"),
                   row_order=NULL,
                   column_order=NULL,
                   column_names_rot = 45,
                   row_names_side = 'right',
                   show_column_names = F,
                   show_row_names = T,  
                   row_names_gp = gpar(col = 'grey30',fontsize = 7),
                   row_names_max_width = unit(12, "cm"),
                   column_labels = colnames(plotdata),
                   column_split =  unique(cc$group) ,
                   row_gap = unit(c(2), "mm"),
                   column_gap = unit(c(2), "mm"),
                   column_title = NULL,
                   column_title_gp = gpar(fill = "red", col = "white", border = "blue",fontface = "bold",fontsize = 9),
                   row_title_gp = gpar(fontsize=12),
                   show_heatmap_legend = T,
                   heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), title_gp = gpar(fontsize = 10, fontface = "bold"))) 





####################################################
######################   TF      ###################
####################################################

exp = mod_data%>%select(-c('group'))
keygene = read.table('TF.txt')
exp_key = as.data.frame(exp) %>% rownames_to_column() %>%
  filter( rowname %in% keygene$V1) %>% 
  column_to_rownames('rowname')


STEMI <- rownames(group_list[group_list$group=='STEMI',]);CCS <- rownames(group_list[group_list$group=='CCS',]);mydata <- exp_key[,c(STEMI,CCS)];group <- factor( c( rep('STEMI',111),rep('CCS',46) ))
design <- model.matrix(~0+group)
colnames(design) <- levels(group)
rownames(design) <- colnames(mydata)
contrst <- makeContrasts('STEMI-CCS',levels = design)

fit <- lmFit(mydata,design = design)
fit2 <- contrasts.fit(fit,contrasts = contrst)
fit2 <- eBayes(fit2)
alldiff <- topTable(fit2,coef = 1,Inf)
sigdiff <- subset(alldiff,
                  abs(logFC) > 0.5 & 
                    adj.P.Val < 0.05)
sig_gene <- rownames(sigdiff)
write.table(sig_gene, file = 'sig_TF.txt')
exp_key = exp_key[sig_gene ,]
write.csv(exp_key,'TF_sig_exp.csv')


#########==========heatmap
hmExp =  read.csv("TF_sig_exp.csv",row.names = 1)
plotdata <- as.matrix(hmExp)

plotdata_mean = cbind(t(standarize.fun(plotdata)),data.frame(group = cc$group)) %>%
                 group_by(group) %>% summarise_all(mean) %>% column_to_rownames('group')  %>% t() %>% as.data.frame()

Left = rowAnnotation( Level = rep('TF', length( row.names(plotdata)))  ,
                      annotation_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"),ncol=1), border = T,
                      col=list(Level =   c( 'TF' =  pal_npg("nrc", alpha = 0.6)(9)[2])),
                      show_annotation_name = F,annotation_name_gp = gpar(fontsize = 10),
                      simple_anno_size_adjust = T,show_legend = T,annotation_width = unit(0.2, "cm")
)

hm4_mean = Heatmap(as.matrix(plotdata_mean),
                   name='TF',
                   left_annotation = Left,
                   cluster_rows = T,
                   cluster_columns = F,
                   color_space = "RGB",
                   col = colorRamp2(seq(-1.19, 1.19, length = 3), c("#84D2EE", 'white','#fec89a')),
                   border = T,
                   border_gp = gpar(col = "black",lty = 1,lwd=1.5, linejoin = 'round' ),
                   rect_gp = gpar(col = "white",lty = 1,lwd=2, linejoin = 'round' ),
                   row_dend_side = c("left"),
                   row_dend_width = unit(1, "mm"),
                   row_order=NULL,
                   column_order=NULL,
                   column_names_rot = 45,
                   row_names_side = 'right',
                   show_column_names = F,
                   show_row_names = T,  
                   row_names_gp = gpar(col = 'grey30',fontsize = 7),
                   row_names_max_width = unit(12, "cm"),
                   column_labels = colnames(plotdata),
                   column_split =  unique(cc$group) ,
                   row_gap = unit(c(2), "mm"),
                   column_gap = unit(c(2), "mm"),
                   column_title = NULL,
                   column_title_gp = gpar(fill = "red", col = "white", border = "blue",fontface = "bold",fontsize = 9),
                   row_title_gp = gpar(fontsize=12),
                   show_heatmap_legend = T,
                   heatmap_legend_param=list(labels_gp = gpar(fontsize = 10), 
                                             title_gp = gpar(fontsize = 10, fontface = "bold"))) 




####################################################
#################   Pathway      ###################
####################################################

exp = mod_data%>%select(-c('group'))
set.seed(1211)
BIOCARTAmarkers <- getGmt('c2.cp.biocarta.v7.5.1.symbols.gmt', geneIdType=SymbolIdentifier())
original <- gsva(as.matrix(exp), BIOCARTAmarkers )
gsva_es = original
row.names(gsva_es) <- str_replace(row.names(gsva_es), "BIOCARTA_", "")
row.names(gsva_es) <- gsub( "PATHWAY", "pathway",row.names(gsva_es))
row.names(gsva_es) <- gsub( "_", " ",row.names(gsva_es))
write.csv(gsva_es, "gsva_output_biocarta.csv", quote = F)

design <- model.matrix(~ 0 + group_list$group)
colnames(design) <- levels(group_list$group)
rownames(design) <- colnames(gsva_es)
fit <- lmFit(gsva_es, design)
fit2 <- eBayes(fit)
alldiff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
write.csv(alldiff, "gsva_limma_biocarta.csv", quote = F)


gsva_es <- read.csv("gsva_output_biocarta.csv", row.names = 1,check.names = F)
alldiff = read.csv( 'gsva_limma_biocarta.csv',row.names = 1 )
sigdiff <- subset(alldiff, abs(t) >= 5 &   adj.P.Val < 0.05)
sig_gene <- rownames(sigdiff) 
write.table(sig_gene, file = 'sig_biocarta.txt')   
exp_key = gsva_es[sig_gene ,] 
write.csv(exp_key,'biocarta_key.csv')


#########==========heatmap
hmExp = read.csv( "biocarta_key.csv",row.names = 1  )
plotdata <- as.matrix(hmExp)

plotdata_mean = cbind(t(standarize.fun(plotdata)),data.frame(group = cc$group)) %>%
  group_by(group) %>% summarise_all(mean) %>% column_to_rownames('group')  %>% t() %>% as.data.frame()

Left = rowAnnotation( Level = rep('PATHWAY', length( row.names(plotdata)))  ,
                      annotation_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"),ncol=1),border = T,
                      col=list(Level = c( 'PATHWAY' =  pal_npg("nrc", alpha = 0.6)(9)[5])),
                      show_annotation_name = F,annotation_name_gp = gpar(fontsize = 10),
                      simple_anno_size_adjust = T,  show_legend = T,annotation_width = unit(0.2, "cm")
)


hm5_mean = Heatmap(as.matrix(plotdata_mean),
                   name=' PATHWAY ',
                   left_annotation = Left,
                   cluster_rows = T,
                   cluster_columns = F,
                   color_space = "RGB",
                   col = colorRamp2(seq(-1.5, 1.5, length = 3), c("#84D2EE", 'white','#FFB562')),
                   border = T,
                   border_gp = gpar(col = "black",lty = 1,lwd=1.5, linejoin = 'round' ),
                   rect_gp = gpar(col = "white",lty = 1,lwd=2, linejoin = 'round' ),
                   row_dend_side = c("left"),
                   row_dend_width = unit(1, "mm"),
                   row_order=NULL,
                   column_order=NULL,
                   column_names_rot = 45,
                   row_names_side = 'right',
                   show_column_names = F,
                   show_row_names = T,  
                   row_names_gp = gpar(col = 'grey30',fontsize = 7),
                   row_names_max_width = unit(12, "cm"),
                   column_labels = colnames(plotdata),
                   column_split =  unique(cc$group) ,
                   row_gap = unit(c(2), "mm"),
                   column_gap = unit(c(2), "mm"),
                   column_title = NULL,
                   column_title_gp = gpar(fill = "red", col = "white", border = "blue",fontface = "bold",fontsize = 9),
                   row_title_gp = gpar(fontsize=12),
                   show_heatmap_legend = T,
                   heatmap_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"))) 



####################################################
#################   REACTOME      ##################
####################################################

exp = mod_data%>%select(-c('group'))

set.seed(1211)
reactomemarkers <- getGmt('c2.cp.reactome.v7.5.1.symbols.gmt', geneIdType=SymbolIdentifier())
original <- gsva(as.matrix(exp), reactomemarkers )
gsva_es = original
row.names(gsva_es) <- str_replace(row.names(gsva_es), "REACTOME_", "")
row.names(gsva_es) <- gsub( "_", " ",row.names(gsva_es))
write.csv(gsva_es, "gsva_output_reactome.csv", quote = F)

design <- model.matrix(~ 0 + group_list$group)
colnames(design) <- levels(group_list$group)
rownames(design) <- colnames(gsva_es)
fit <- lmFit(gsva_es, design)
fit2 <- eBayes(fit)
alldiff <- topTable(fit2, coef = 1, n = Inf, adjust.method = "BH", sort.by = "P")
write.csv(alldiff, "gsva_limma_reactome.csv", quote = F)

alldiff = read.csv("gsva_limma_reactome.csv", row.names = 1)
sigdiff <- subset(alldiff, abs(t) >= 5 &  
                    adj.P.Val < 0.05)
sig_gene <- rownames(sigdiff) 
write.table(sig_gene, file = 'sig_reactome.txt')   
exp_key = gsva_es[sig_gene ,] 
write.csv(exp_key,'reactome_key.csv')


#########==========heatmap
hmExp = read.csv( "reactome_key.csv",row.names = 1  )
plotdata <- as.matrix(hmExp)

plotdata_mean = cbind(t(standarize.fun(plotdata)),data.frame(group = cc$group)) %>%
      group_by(group) %>% summarise_all(mean) %>% column_to_rownames('group')  %>% t() %>% as.data.frame()

Left = rowAnnotation( Level = rep('REACTOME', length( row.names(plotdata)))  ,
                      annotation_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"),
                                                   ncol=1),border = T,
                      col=list(Level = c('REACTOME' =  pal_npg("nrc", alpha = 0.6)(9)[6])),
                      show_annotation_name = F,annotation_name_gp = gpar(fontsize = 10),
                      simple_anno_size_adjust = T, show_legend = T,annotation_width = unit(0.2, "cm")
)

hm6_mean = Heatmap(as.matrix(plotdata_mean),
                   name=' REACTOME ',
                   left_annotation = Left,
                   cluster_rows = T,
                   cluster_columns = F,
                   color_space = "RGB",
                   col = colorRamp2(seq(-1.19, 1.19, length = 3), c("#84D2EE", 'white','#ffadad')),
                   border = T,
                   border_gp = gpar(col = "black",lty = 1,lwd=1.5, linejoin = 'round' ),
                   rect_gp = gpar(col = "white",lty = 1,lwd=2, linejoin = 'round' ),
                   row_dend_side = c("left"),
                   row_dend_width = unit(1, "mm"),
                   row_order=NULL,
                   column_order=NULL,
                   column_names_rot = 45,
                   row_names_side = 'right',
                   show_column_names = F,
                   show_row_names = T,  
                   row_names_gp = gpar(col = 'grey30',fontsize = 7),
                   row_names_max_width = unit(12, "cm"),
                   column_labels = colnames(plotdata),
                   column_split =  unique(cc$group) ,
                   row_gap = unit(c(2), "mm"),
                   column_gap = unit(c(2), "mm"),
                   column_title = NULL,
                   column_title_gp = gpar(fill = "red", col = "white", border = "blue",fontface = "bold",fontsize = 9),
                   row_title_gp = gpar(fontsize=12),
                   show_heatmap_legend = T,
                   heatmap_legend_param=list(labels_gp = gpar(fontsize = 10),title_gp = gpar(fontsize = 10, fontface = "bold"))) 



####################################################
#############   Complexheatmap      ################
####################################################

draw(hm2_mean %v% hm4_mean %v% hm3_mean %v% hm5_mean %v% hm1_mean %v% hm6_mean, 
     heatmap_legend_side = "bottom", annotation_legend_side = "bottom")


####################################################
###############   Interaction      #################
####################################################


#########==========Cor HALLMARK

TF <- mod_data %>% select(-c("group"))
uniSigExp <- read.table("APVS gene.txt") %>% .[, 1]
gene1 <- uniSigExp$V1
TF <- TF[gene1, ]
gsva_es <- read.csv("gsva_output_hall.csv", row.names = 1)

diff <- read.table("sig_HALLMARK.txt", row.names = 1)
diffName <- diff$x
gsva_es <- gsva_es[diffName, ]
immuneGene <- gsva_es

sameSample <- intersect(colnames(TF), colnames(immuneGene))
TF1 <- TF[, sameSample]
immuneGene1 <- immuneGene[, sameSample]

corFilter <- 0.6
pvalueFilter <- 0.0001
outTab <- data.frame()

for (i in row.names(TF1)) {
  for (j in row.names(immuneGene1)) {
    x <- as.numeric(TF1[i, ])
    y <- as.numeric(immuneGene1[j, ])
    corT <- cor.test(x, y)
    cor <- corT$estimate
    pvalue <- corT$p.value
    if ((cor > corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "postive"))
    }
    if ((cor < -corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "negative"))
    }
  }
}
outTab <- outTab %>% filter(cor != 1)
write.table(file = "corResult_HALL.txt", outTab, sep = "\t", quote = F, row.names = F)



#########==========Cor TF

TF <- mod_data %>% select(-c("group"))
uniSigExp <- read.table("APVS gene.txt") %>% .[, 1]
gene1 <- uniSigExp$V1
TF <- TF[gene1, ]

immuneGene <- read.csv("TF_sig_exp.csv", row.names = 1)
sameSample <- intersect(colnames(TF), colnames(immuneGene))
TF1 <- TF[, sameSample]
immuneGene1 <- immuneGene[, sameSample]

corFilter <- 0.6
pvalueFilter <- 0.0001
outTab <- data.frame()
for (i in row.names(TF1)) {
  for (j in row.names(immuneGene1)) {
    x <- as.numeric(TF1[i, ])
    y <- as.numeric(immuneGene1[j, ])
    corT <- cor.test(x, y)
    cor <- corT$estimate
    pvalue <- corT$p.value
    if ((cor > corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "postive"))
    }
    if ((cor < -corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "negative"))
    }
  }
}
outTab <- outTab %>% filter(cor != 1)
write.table(file = "corResult_TF.txt", outTab, sep = "\t", quote = F, row.names = F)



#########==========Cor cell

TF <- mod_data %>% select(-c("group"))
uniSigExp <- read.table("APVS gene.txt") %>% .[, 1]
gene1 <- uniSigExp$V1
TF <- TF[gene1, ]

exp <- read.csv("xCell.csv", header = T, row.names = 1)
sig_gene <- read.table("sig_xcell.txt")[, 1]
immuneGene <- exp[sig_gene, ]

sameSample <- intersect(colnames(TF), colnames(immuneGene))
TF1 <- TF[, sameSample]
immuneGene1 <- immuneGene[, sameSample]

corFilter <- 0.6
pvalueFilter <- 0.0001
outTab <- data.frame()
for (i in row.names(TF1)) {
  for (j in row.names(immuneGene1)) {
    x <- as.numeric(TF1[i, ])
    y <- as.numeric(immuneGene1[j, ])
    corT <- cor.test(x, y)
    cor <- corT$estimate
    pvalue <- corT$p.value
    if ((cor > corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "postive"))
    }
    if ((cor < -corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "negative"))
    }
  }
}
outTab <- outTab %>% filter(cor != 1)
write.table(file = "corResult_xCell.txt", outTab, sep = "\t", quote = F, row.names = F)



#########==========Cor PATHWAY

TF <- mod_data %>% select(-c("group"))
uniSigExp <- read.table("APVS gene.txt") %>% .[, 1]
gene1 <- uniSigExp$V1
TF <- TF[gene1, ]

gsva_es <- read.csv("gsva_output_biocarta.csv", row.names = 1)
diff <- read.table("sig_biocarta.txt", row.names = 1)
diffName <- diff$x
gsva_es <- gsva_es[diffName, ]
immuneGene <- gsva_es
sameSample <- intersect(colnames(TF), colnames(immuneGene))
TF1 <- TF[, sameSample]
immuneGene1 <- immuneGene[, sameSample]

corFilter <- 0.6
pvalueFilter <- 0.0001
outTab <- data.frame()
for (i in row.names(TF1)) {
  for (j in row.names(immuneGene1)) {
    x <- as.numeric(TF1[i, ])
    y <- as.numeric(immuneGene1[j, ])
    corT <- cor.test(x, y)
    cor <- corT$estimate
    pvalue <- corT$p.value
    if ((cor > corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "postive"))
    }
    if ((cor < -corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "negative"))
    }
  }
}
outTab <- outTab %>% filter(cor != 1)
write.table(file = "corResult_biocarta.txt", outTab, sep = "\t", quote = F, row.names = F)



#########==========Cor reactome

TF <- mod_data %>% select(-c("group"))
uniSigExp <- read.table("APVS gene.txt") %>% .[, 1]
gene1 <- uniSigExp$V1
TF <- TF[gene1, ]

gsva_es <- read.csv("gsva_output_reactome.csv", row.names = 1)
diff <- read.table("sig_reactome.txt", row.names = 1)
diffName <- diff$x
gsva_es <- gsva_es[diffName, ]
immuneGene <- gsva_es
sameSample <- intersect(colnames(TF), colnames(immuneGene))
TF1 <- TF[, sameSample]
immuneGene1 <- immuneGene[, sameSample]

corFilter <- 0.6
pvalueFilter <- 0.0001
outTab <- data.frame()
for (i in row.names(TF1)) {
  for (j in row.names(immuneGene1)) {
    x <- as.numeric(TF1[i, ])
    y <- as.numeric(immuneGene1[j, ])
    corT <- cor.test(x, y)
    cor <- corT$estimate
    pvalue <- corT$p.value
    if ((cor > corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "postive"))
    }
    if ((cor < -corFilter) & (pvalue < pvalueFilter)) {
      outTab <- rbind(outTab, cbind(TF = i, immuneGene = j, cor, pvalue, Regulation = "negative"))
    }
  }
}
outTab <- outTab %>% filter(cor != 1)
write.table(file = "corResult_reactome.txt", outTab, sep = "\t", quote = F, row.names = F)




######### ==========Cor ALL
cor1 <- read.table("corResult_TF.txt", header = T, sep = "\t", check.names = F)
cor2 <- read.table("corResult_HALL.txt", header = T, sep = "\t", check.names = F)
cor3 <- read.table("corResult_xCell.txt", header = T, sep = "\t", check.names = F)
cor4 <- read.table("corResult_biocarta.txt", header = T, sep = "\t", check.names = F)
cor5 <- read.table("corResult_reactome.txt", header = T, sep = "\t", check.names = F)

lnc <- intersect(as.vector(cor1[, 1]), as.vector(cor2[, 1]))
lnc <- intersect(lnc, as.vector(cor3[, 1]))
lnc <- intersect(lnc, as.vector(cor4[, 1]))
lnc <- intersect(lnc, as.vector(cor5[, 1]))
write.table(file = "final_cor_sig.txt", lnc, sep = "\t", quote = F, row.names = F)

cor6 <- cor1[which(cor1[, 1] %in% lnc), ]
cor7 <- cor2[which(cor2[, 1] %in% lnc), ]
cor8 <- cor3[which(cor3[, 1] %in% lnc), ]
cor9 <- cor4[which(cor4[, 1] %in% lnc), ]
cor10 <- cor5[which(cor5[, 1] %in% lnc), ]

cor_all=rbind(cor6,cor7,cor8,cor9,cor10)
write.table(file="cor_all.txt",cor_all,sep="\t",quote=F,row.names=F) 
cor_all <- read.table(file="cor_all.txt",sep="\t")
write.csv(cor_all,file = "cor_all.csv")


listgene = c(levels(factor(cor6$TF)),levels(factor(cor6$immuneGene)))
listgene = unique(listgene)  
write.table(file="finallistgene.txt",listgene,sep="\t",quote=F,row.names=F) 
listhallmark = c(levels(factor(cor7$immuneGene)))
write.table(file="finallisthallmark.txt",listhallmark,sep="\t",quote=F,row.names=F) 
listssgsea = c(levels(factor(cor8$immuneGene)))
write.table(file="finallistxcell.txt",listssgsea,sep="\t",quote=F,row.names=F) 
listpathway = c(levels(factor(cor9$immuneGene)))
write.table(file="finallistbiocarta.txt",listpathway,sep="\t",quote=F,row.names=F) 
listreactome = c(levels(factor(cor10$immuneGene)))
write.table(file="finallistreactome.txt",listreactome,sep="\t",quote=F,row.names=F) 

nodeType6 = data.frame(symbol = levels(factor(cor6$immuneGene)), type = 'TF')  
nodeType6_2 = data.frame(symbol = levels(factor(cor6$TF)), type = 'APVS') 
nodeType7 = data.frame(symbol = levels(factor(cor7$immuneGene)), type = 'Hallmark')
nodeType8 = data.frame(symbol = levels(factor(cor8$immuneGene)), type = 'Cell Landscape')
nodeType9 = data.frame(symbol = levels(factor(cor9$immuneGene)), type = 'Pathway')
nodeType10 = data.frame(symbol = levels(factor(cor10$immuneGene)), type = 'Reactome')

nodeType = rbind(nodeType6,nodeType6_2,nodeType7,nodeType8, nodeType9 , nodeType10 )
write.table(file="nodeType.txt",nodeType,sep="\t",quote=F,row.names=F) 





gsva_es <- read.csv("gsva_output_hall.csv", row.names = 1)
TIMER2 <- read.csv("xCell.csv", header = T, row.names = 1)
rt <- mod_data %>% select(-c("group"))
rt <- as.matrix(rt)
rownames(rt) <- rt[, 1]
exp <- rt[, 2:ncol(rt)]
dimnames <- list(rownames(exp), colnames(exp))
immuneGene <- matrix(as.numeric(as.matrix(exp)), nrow = nrow(exp), dimnames = dimnames)

gsva_biocarta <- read.csv("gsva_output_biocarta.csv", row.names = 1)
gsva_reactome <- read.csv("gsva_output_reactome.csv", row.names = 1)

listgene <- read.table("finallistgene.txt", header = T)
listpathway <- read.table("finallistbiocarta.txt", header = T, sep = "\t")
listssgsea <- read.table("finallistxcell.txt", header = T, sep = "\t")
listhallmark <- read.table("finallisthallmark.txt", header = T, sep = "\t")
listreactome <- read.table("finallistreactome.txt", header = T, sep = "\t")

listgene <- as.vector(listgene[, 1])
listpathway <- as.vector(listpathway[, 1])
listssgsea <- as.vector(listssgsea[, 1])
listhallmark <- as.vector(listhallmark[, 1])
listreactome <- as.vector(listreactome[, 1])

finalgene <- immuneGene[listgene, ]
finalhallmark <- gsva_es[listhallmark, ]
finalssgeas <- TIMER2[listssgsea, ]
finalpathway <- gsva_biocarta[listpathway, ]
finalreactome <- gsva_reactome[listreactome, ]

finalall <- rbind(finalgene, finalpathway, finalssgeas, finalhallmark, finalreactome)

ggcorrplot(
  corr = cor(t(finalall)),
  type = "lower",
  ggtheme = ggplot2::theme_classic(base_rect_size = 1),
  colors = c("#4D90D4", "white", "#E9536B"),
  lab = F, lab_size = 0.9, tl.cex = 1
) + labs(y = "", x = "")

corfinal <- cor(t(finalall))
write.csv(corfinal, 'corfinal.csv')







