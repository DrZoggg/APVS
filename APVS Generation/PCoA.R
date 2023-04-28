
if (T) {
library(vegan)        
library(ape)          
library(RColorBrewer) 
library(usedist)      
library(ggsci)
library("scales")
}

APVS <- read.table("APVS gene.txt") %>% .[, 1]
expr1 <- read.csv("GSE59867_expr.csv", row.names = 1, header = T)
expr1 <- expr1[APVS, ]
sinfo1 <- read.csv("GSE59867_clin.csv", row.names = 1)
expr2 <- read.csv("GSE62646_expr.csv", row.names = 1, header = T)
expr2 <- expr2[APVS, ]
sinfo2 <- read.csv("GSE62646_clin.csv", row.names = 1)

expr <- cbind(expr1, expr2)
sinfo1 <- sinfo1[, 3, drop = F]
sinfo1$Cohort <- rep("GSE59867", 436)
sinfo2 <- sinfo2[, 3, drop = F]
sinfo2$Cohort <- rep("GSE62646", 98)
sinfo <- rbind(sinfo1, sinfo2)


####################################################
####################   PCoA      ###################
####################################################

brca.dist <- vegdist(t(expr), method = "euclidean")
pc <- cmdscale(brca.dist, eig = FALSE)
eig <- pcoa(brca.dist)
pct1 <- round(eig$values$Relative_eig[1] / sum(eig$values$Relative_eig), digits = 3) * 100
pct2 <- round(eig$values$Relative_eig[2] / sum(eig$values$Relative_eig), digits = 3) * 100
xlab.text <- paste("PCoA1 [", pct1, "%]", sep = "")
ylab.text <- paste("PCoA2 [", pct2, "%]", sep = "")


###### ===plot1
plotinfo <- cbind.data.frame(
  x = pc[, 1],
  y = pc[, 2],
  group = sinfo[rownames(pc), "group"],
  Cohort = sinfo[rownames(pc), "Cohort"],
  Class = paste(sinfo[rownames(pc), "group"],
    sinfo[rownames(pc), "Cohort"],
    sep = "_"
  )
)

plotinfo$group <- factor(plotinfo$group,
  levels = c(
    "Stable CCS", "1st day of STEMI", "4-6 days of STEMI",
    "1 month after STEMI", "6 months after STEMI"
  )
)

mypal <- pal_uchicago(alpha = .8)(10)
plotinfo$color <- mypal[plotinfo$group]

par(bty = "o", mgp = c(1.9, .33, 0), mar = c(3.1, 3.1, 2.1, 2.1) + .1, las = 1, tcl = -.25)
plot(plotinfo$x,
  plotinfo$y,
  pch = 19,
  col = plotinfo$color,
  xlab = xlab.text,
  ylab = ylab.text
)


######===plot2
plotinfo2 <- NULL
for (i in unique(plotinfo$Class)) {
  tmp <- plotinfo[plotinfo$Class == i, ]
  avgx <- mean(tmp$x)
  avgy <- mean(tmp$y)
  sdx <- sd(tmp$x)
  sdy <- sd(tmp$y)

  plotinfo2 <- rbind.data.frame(plotinfo2,
    data.frame(
      Class = i,
      color = unique(tmp$color),
      shape = ifelse(unique(tmp$Cohort) == "GSE62646",
        "closed", "opened"
      ),
      label = switch(i,
        "1st day of STEMI_GSE59867" = "STEMI1D",
        "4-6 days of STEMI_GSE59867" = "4-6D",
        "1 month after STEMI_GSE59867" = "1Mth",
        "6 months after STEMI_GSE59867" = "6Mth",
        "Stable CCS_GSE59867" = "CCS",
        "1st day of STEMI_GSE62646" = "STEMI1D",
        "4-6 days of STEMI_GSE62646" = "4-6D",
        "6 months after STEMI_GSE62646" = "6Mht",
        "Stable CCS_GSE62646" = "CCS"
      ),
      avgx = avgx,
      avgy = avgy,
      sdx = sdx,
      sdy = sdy,
      stringsAsFactors = F
    ),
    stringsAsFactors = F
  )
}


par(bty="o", mgp = c(1.9,.33,0), mar=c(3.1,3.1,2.1,2.1)+.1, las=1, tcl=-.25)
plot(NULL,NULL,
     xlab = xlab.text,
     ylab = ylab.text,
     xlim = c(-4,2.3),
     axes = T,
     mar=c(4, 4, 2, 2)+.1,
     mgp=c(1.5, .2, 5),
     ylim = c(-0.5,0.55))

for (i in 1:nrow(plotinfo2)) {
  tmp <- plotinfo2[i,]
  
  lines(x = c(tmp$avgx - tmp$sdx,
              tmp$avgx + tmp$sdx),
        y = c(tmp$avgy, tmp$avgy),
        lty = ifelse(tmp$shape == "closed",1,2), 
        col = tmp$color,
        lwd = 2)
  
  lines(x = c(tmp$avgx - tmp$sdx,
              tmp$avgx - tmp$sdx),
        y = c(tmp$avgy, 
              tmp$avgy),
        col = tmp$color,
        lwd = 2)
  lines(x = c(tmp$avgx + tmp$sdx,
              tmp$avgx + tmp$sdx),
        y = c(tmp$avgy,
              tmp$avgy),
        col = tmp$color,
        lwd = 2)
  
  lines(x = c(tmp$avgx, tmp$avgx),
        y = c(tmp$avgy - tmp$sdy,
              tmp$avgy + tmp$sdy),
        lty = ifelse(tmp$shape == "closed",1,2),
        col = tmp$color,
        lwd = 2)
  
  lines(x = c(tmp$avgx,
              tmp$avgx),
        y = c(tmp$avgy + tmp$sdy,
              tmp$avgy + tmp$sdy),
        col = tmp$color,
        lwd = 2)
  lines(x = c(tmp$avgx,
              tmp$avgx),
        y = c(tmp$avgy - tmp$sdy,
              tmp$avgy - tmp$sdy),
        col = tmp$color,
        lwd = 2)
}

points(plotinfo2$avgx,
       plotinfo2$avgy,
       pch = ifelse(plotinfo2$shape == "closed",19,21), 
       bg = ifelse(plotinfo2$shape == "closed",plotinfo2$color,"white"),
       col = plotinfo2$color, 
       lwd = 2, 
       cex = 5.5) 

text(plotinfo2$avgx,
     plotinfo2$avgy,
     plotinfo2$label,
     col = ifelse(plotinfo2$shape == "closed","white",plotinfo2$color), 
     cex = 0.8)




