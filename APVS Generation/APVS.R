

if (T) {
  library(tidyverse)
  library(dplyr)
  library(data.table)
  library(DALEX)
  library(ggplot2)
  library(randomForest)
  library(kernlab)
  library(caret)
  library(pROC)
  library(e1071)
  library(glmnet)
  library(randomForest)
  library(mboost)
  library(kknn)
  library(pls)
  library(rpart)
  library(neuralnet)
  library(patchwork)
  library(reshape2)
  library(scales)
  library(ggpubr)
  library(pheatmap)
  library(ComplexHeatmap)
  library(data.table)
  library(Matrix)
  library(viridis)
  library(harmony)
  library(RColorBrewer)
  library(tidyverse)
  library(dplyr)
  library(ggthemr)
  library(gghalves)
  library(scales)
  library(ggsci)
}



DCPGs <- fread("Dysregulated gene co-expression pattern.csv", header = T) %>% pull(DCPGs)
mod_data <- fread("expr.csv", header = T) %>%
  filter(V1 %in% DCPGs) %>%
  column_to_rownames("V1") %>%
  t() %>%
  merge(., read.csv("clin.csv", row.names = 1, header = T) %>% select(group), by.x = 0, by.y = 0) %>%
  column_to_rownames("Row.names")

mod <- list()
if (T) {
  set.seed(1000)
  index <- createDataPartition(mod_data$Type, p = 0.7, list = F)
  train_data <- mod_data[index, ]
  test_data <- mod_data[-index, ]

  set.seed(1000)
  control <- trainControl(method = "repeatedcv", number = 10, savePredictions = TRUE, repeats = 100)

  mod_rf <- train(Type ~ .,
    data = train_data,
    ntree = 3000,
    method = "rf",
    trControl = control,
    tuneLength = 10
  );mod[[mod_rf]] <- mod_rf

  mod_treebag <- train(Type ~.,
                      data= train_data,
                      method = "treebag",
                      trControl = control,
                      tuneLength = 10
  );mod[[mod_treebag]] <- mod_treebag

  mod_pls <- train(Type ~ .,
    data = train_data,
    method = "pls",
    prob.model = TRUE,
    trControl = control,
    tuneLength = 10
  );mod[[mod_pls]] <- mod_pls

  mod_nnet <- train(Type ~ .,
    data = train_data,
    method = "nnet",
    prob.model = TRUE,
    preProc = c("range"), rangeBounds = c(0, 1),
    trControl = control,
    tuneLength = 10
  );mod[[mod_nnet]] <- mod_nnet

  mod_glmnet <- train(Type ~ .,
    data = train_data,
    method = "glmnet",
    trControl = control,
    tuneLength = 10
  );mod[[mod_glmnet]] <- mod_glmnet

  mod_bayes <- train(Type ~ .,
    data = train_data,
    method = "naive_bayes",
    trControl = control,
    tuneLength = 10
  );mod[[mod_bayes]] <- mod_bayes

  mod_knn <- train(Type ~ .,
    data = train_data,
    method = "knn",
    trControl = control,
    tuneLength = 10
  );mod[[mod_knn]] <- mod_knn

  mod_rpart <- train(Type ~ .,
    data = train_data,
    method = "rpart",
    trControl = control,
    tuneLength = 10
  );mod[[mod_rpart]] <- mod_rpart

  mod_glmboost <- train(Type ~ .,
    data = train_data,
    method = "glmboost",
    trControl = control,
    tuneLength = 10
  );mod[[mod_glmboost]] <- mod_glmboost
}


####################################################
################   Evaluation       ################
####################################################


p_fun=function(object, newdata){
  predict(object, newdata=newdata, 
          type="prob"
  )[,2]
}
yTest = ifelse(test_data$Type=="CAD", 0, 1)

performance <- data.frame()
mp_all <- list()
for (i in names(mod)) {
  explainer <- explain(mod[[i]],
    data = test_data, y = yTest,
    predict_function = p_fun,
    verbose = T
  )
  mp <- model_performance(explainer)
  mp_all[[i]] <- mp
  performance <- cbind(performance, unlist(mp$measures) %>% data.frame())
}
save(performance, file = "performance.rda")



combined_df <- reduce(mp_all, rbind)
residual_dataframe <- data.frame(label = c(1:46))
label <- c("rf", "knn", "nnet", "bayes", "glmnet", "rpart", "glmboost", "pls", "treebag")
for (i in 1:length(label)) {
  residual_dataframe <- cbind(
    residual_dataframe,
    data.frame(name = abs(as.numeric(combined_df[i, 1][[1]]$diff)))
  )
  colnames(residual_dataframe)[i + 1] <- label[i]
}
residual_dataframe <- residual_dataframe %>%
  select(-1) %>%
  as.data.frame() %>%
  gather(key = model_type, value = residual)

RMSE_dataframe <- data.frame(
  RMSE = residual_dataframe %>%
    group_by(model_type) %>%
    summarise(RMSE = sqrt(mean(residual^2))) %>%
    arrange(desc(RMSE)) %>%
    pull(RMSE),
  model_type = residual_dataframe %>%
    group_by(model_type) %>%
    summarise(RMSE = sqrt(mean(residual^2))) %>%
    arrange(desc(RMSE)) %>%
    pull(model_type) %>%
    as.character()
)
RMSE_residual <- merge(RMSE_dataframe, residual_dataframe, by = "model_type")
RMSE_residual$model_type <- factor(RMSE_residual$model_type, levels = RMSE_dataframe$model_type)
save(RMSE_residual, file = "RMSE_residual.rda")



####################################################
##############   Evaluation plot      ##############
####################################################

load("RMSE_residual.rda")
load("performance.rda")

performance_ggdata <- performance %>%
  as.data.frame() %>%
  gather(key = model_type, value = value, -X)
colnames(performance_ggdata) <- c("index", "model", "value")

ggthemr("fresh", line_weight = .3)
p1 <- ggplot(
  performance_ggdata,
  aes(x = value, y = model, colour = model)
) +
  facet_wrap(~index, ncol = 9) +
  scale_colour_uchicago(alpha = .6) +
  geom_point(size = 4, stroke = 1.5, shape = 21) +
  theme_bw(base_rect_size = 1.2) +
  theme(
    panel.grid.minor.y = element_blank(),
    panel.grid.major.y = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold", angle = 90),
    axis.title = element_blank(),
    axis.ticks.y = element_blank(),
    axis.line.y = element_blank(),
    axis.text.y = element_text(size = 10, face = "bold", colour = "#41140C"),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    strip.background = element_rect(),
    strip.text = element_blank()
  ) +
  guides(size = "none")


p2 <- ggplot(RMSE_residual) +
  geom_boxplot(aes(x = model_type, y = residual, fill = model_type),
    width = 0.15,
    alpha = 0.6,
    position = position_dodge(0.2),
    outlier.shape = NA
  ) +
  coord_flip() +
  scale_fill_uchicago(alpha = 1) +
  geom_half_dotplot(aes(x = model_type, y = residual, fill = model_type),
    alpha = 0.4,
    method = "histodot", stackdir = "down",
    dotsize = 1.2,
    position = position_nudge(x = -0.15, y = 0),
    binwidth = 0.03,
    colour = NA
  ) +
  geom_point(aes(x = model_type, y = RMSE),
    size = 4, shape = 18,
    color = alpha("#f95738", 0.4)
  ) +
  theme_classic(base_rect_size = 1.2) +
  theme(
    axis.title.y = element_blank(),
    axis.text.x = element_text(size = 9, face = "bold", angle = 90),
    axis.title.x = element_blank(),
    axis.line.y = element_blank(),
    axis.line.x = element_blank(),
    axis.text.y = element_blank(),
    axis.ticks.y = element_blank(),
    legend.title = element_text(face = "bold", size = 10),
    axis.ticks.x = element_blank(),
    legend.position = "none",
    legend.key.size = unit(6, "mm")
  )


ggarrange(p1, p2,
  nrow = 1, ncol = 2,
  align = "v", widths = c(9, 3), common.legend = F
)




####################################################
##############   APVS generation      ##############
####################################################

DCPGs <- fread("Dysregulated gene co-expression pattern.csv", header = T) %>% pull(DCPGs)
mod_data <- fread("expr.csv", header = T) %>%
  filter(V1 %in% DCPGs) %>%
  column_to_rownames("V1") %>%
  t() %>%
  merge(., read.csv("clin.csv", row.names = 1, header = T) %>% select(group), by.x = 0, by.y = 0) %>%
  column_to_rownames("Row.names")

set.seed(1000)
mod <- randomForest(as.factor(group) ~ .,
  data = mod_data,
  ntree = 3000,
  importance = T,
  proximity = T
)
set.seed(1000)
optionTrees <- which.min(rf$err.rate[, 1])
mod <- randomForest(as.factor(group) ~ .,
  data = mod_data,
  ntree = optionTrees,
  importance = T,
  proximity = T
)

## ==importance
raw.imp <- mod$importance
names(raw.imp) <- gsub("_", "-", names(raw.imp))
raw.imp <- data.frame(raw.imp)

rel.imp <- data.frame(raw.imp / max(abs(raw.imp)))
rel.imp <- data.frame(rel.imp)

imp.res <- data.frame(
  gene = row.names(raw.imp),
  raw.importance = raw.imp,
  rel.importance = rel.imp,
  stringsAsFactors = F
)
imp.cutoff <- 0.0
rel.imp.sel <- rel.imp[rel.imp$MeanDecreaseGini > imp.cutoff, , drop = F]
rel.imp.sel <- rel.imp.sel[order(rel.imp.sel$MeanDecreaseGini, decreasing = F), ]


importance <- importance(x = mod)
df <- as.data.frame(importance) %>%
  select(MeanDecreaseGini, MeanDecreaseAccuracy) %>%
  rownames_to_column("id")
my_palette <- colorRampPalette(rev(paletteer_c("grDevices::RedOr", 30)), alpha = .1)(n = 30)


ggplot(df, aes(MeanDecreaseGini, reorder(id, MeanDecreaseGini))) +
  geom_point(aes(x = MeanDecreaseAccuracy),
    shape = 21, size = 2, colour = "white", fill = "white",
    stroke = 0.8, alpha = .6
  ) +
  geom_line(aes(x = MeanDecreaseAccuracy), size = 1.5, group = 1, colour = "#a2d2ff", alpha = .6) +
  geom_bar(
    stat = "identity",
    width = 0.7,
    color = "black",
    aes(fill = MeanDecreaseGini),
    size = 0.2
  ) +
  scale_fill_gradientn(colours = my_palette) +
  theme_bw(base_rect_size = 1.2) +
  ylab(NULL) +
  xlab("") +
  scale_y_discrete(position = "left", ) +
  scale_x_continuous(expand = c(0.01, 0.01)) +
  labs(title = "Variable Importance", color = "FDR") +
  theme(
    axis.text.x = element_text(size = 6, colour = "grey30", face = "bold"),
    axis.text.y = element_text(size = 6, colour = "grey30", face = "bold"),
    axis.title.x = element_blank(),
    axis.title.y = element_text(size = 10, colour = "#41140C", face = "bold"),
    plot.title = element_text(hjust = .5, size = 10, colour = "#41140C", face = "bold"),
    axis.ticks.y = element_blank(),
    panel.grid = element_blank(),
    axis.line.x = element_blank(),
    axis.line.y = element_blank(),
    legend.position = c("none"),
    legend.title = element_text(size = 11, colour = "black"),
    legend.text = element_text(size = 10, colour = "black"),
    legend.background = element_blank(),
    legend.key = element_blank()
  )


Genes <- importance[order(importance[, "MeanDecreaseGini"], decreasing = TRUE), ]
set.seed(1000)
otu_train.cv <- replicate(10, rfcv(mod_data[, -ncol(mod_data)], as.factor(mod_data$group),
  cv.fold = 10, step = 0.5
), simplify = FALSE)
otu_train.cv2 <- data.frame(sapply(otu_train.cv, "[[", "error.cv"))
otu_train.cv2$otus <- rownames(otu_train.cv2)
otu_train.cv2 <- reshape2::melt(otu_train.cv2, id = "otus")
otu_train.cv2$otus <- as.numeric(as.character(otu_train.cv2$otus))

ggplot(otu_train.cv2, aes(otus, value)) +
  geom_smooth(size = 1.8) +
  theme_bw(base_rect_size = 1.2) +
  labs(title = "", x = "Number of Features", y = "Cross-validation error") +
  theme(
    axis.text.x = element_text(size = 10, colour = "black", face = "bold"),
    axis.text.y = element_text(size = 10, colour = "black", face = "bold"),
    axis.title.x = element_text(size = 13, colour = "#41140C", face = "bold"),
    axis.title.y = element_text(size = 13, colour = "#41140C", face = "bold.italic"),
    axis.ticks.y = element_blank(),
    axis.line.y = element_line(size = .8),
    axis.line.x = element_line(size = .8),
    legend.position = c("none"),
    legend.key.size = unit(.2, "cm"),
    legend.background = element_blank(),
    legend.text = element_text(size = 10, colour = "grey40", face = "bold"),
    legend.title = element_text(size = 10, colour = "grey40", face = "bold"),
    plot.title = element_blank()
  )
APVS <- rownames(Genes[1:14, ])



####################################################
##############   APVS classifier      ##############
####################################################

mod_data <- fread("expr.csv", header = T) %>%
  filter(V1 %in% APVS) %>%
  column_to_rownames("V1") %>%
  t() %>%
  merge(., read.csv("clin.csv", row.names = 1, header = T) %>% select(group), by.x = 0, by.y = 0) %>%
  column_to_rownames("Row.names")

set.seed(1000)
index <- createDataPartition(mod_data$group, p = 0.7, list = F)
train_data <- mod_data[index, ]
test_data <- mod_data[-index, ]

DEG_Fun <- function(exprSet,design,contrast.matrix){
  fit <- lmFit(exprSet,design)  ##step1
  fit2 <- contrasts.fit(fit, contrast.matrix) ##step2
  fit2 <- eBayes(fit2) 
  tempOutput = topTable(fit2, coef=1, n=Inf)##step3
  nrDEG = na.omit(tempOutput) 
  return(nrDEG)
}


###=======genesocre

## discovery
group_list <- train_data$group
design <- model.matrix(~ 0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(train_data[, -ncol(train_data)])
contrast.matrix <- makeContrasts(paste("STEMI", "-", "CCS"), levels = design)

data <- train_data[, -ncol(train_data)]
diffRT <- DEG_Fun(data, design, contrast.matrix)

dataUp <- data[diffRT[, "logFC"] > 0, ]
dataDown <- data[diffRT[, "logFC"] < 0, , drop = F]
dataUp2 <- t(apply(dataUp, 1, function(x) ifelse(x > median(x), 1, 0)))
dataDown2 <- t(apply(dataDown, 1, function(x) ifelse(x > median(x), 0, 1)))
outTab <- rbind(dataUp2, dataDown2)

train <- as.data.frame(t(outTab))
train$con <- ifelse(train_data$group == "CCS", 1, 0)
train$treat <- ifelse(train_data$group == "STEMI", 1, 0)

## testing
group_list <- test_data$group
design <- model.matrix(~ 0 + factor(group_list))
colnames(design) <- levels(factor(group_list))
rownames(design) <- colnames(test_data[, -ncol(test_data)])
contrast.matrix <- makeContrasts(paste("STEMI", "-", "CCS"), levels = design)

data <- test_data[, -ncol(test_data)]
diffRT <- DEG_Fun(data, design, contrast.matrix)

dataUp <- data[diffRT[, "logFC"] > 0, ]
dataDown <- data[diffRT[, "logFC"] < 0, , drop = F]
dataUp2 <- t(apply(dataUp, 1, function(x) ifelse(x > median(x), 1, 0)))
dataDown2 <- t(apply(dataDown, 1, function(x) ifelse(x > median(x), 0, 1)))
outTab <- rbind(dataUp2, dataDown2)

test <- as.data.frame(t(outTab))
test$con <- ifelse(test_data$group == "CCS", 1, 0)
test$treat <- ifelse(test_data$group == "STEMI", 1, 0)

save(train, test, file = "genescore.rda")


###=======classifier

load("genescore.rda")
set.seed(1000)
fit <- neuralnet(con + treat ~ ., train,
  hidden = c(10, 6, 3), rep = 10
)
weights <- fit$weights[[1]]
result <- fit$result.matrix

yTest <- test$treat
group <- ifelse(yTest == 1, "STEMI", "CCS")


p_fun <- function(object, newdata) {
  predict(object,
    newdata = newdata,
    type = "prob"
  )[, 2]
}
explainer_fit <- explain(fit,
  label = "nnet",
  data = test[, c(1:14)],
  y = yTest,
  predict_function = p_fun,
  verbose = T
)
mp_fit <- model_performance(explainer_fit)
performance_fit <- data.frame(nnet = unlist(mp_fit$measures))

net.predict <- compute(fit, test[, 1:14], rep = 1)$net.result

net.prediction <- c("CCS", "STEMI")[apply(net.predict, 1, which.max)]
predict.table <- table(group, net.prediction)

classAgreement(predict.table)
confusionMatrix <- confusionMatrix(predict.table)

conAccuracy <- predict.table[1, 1] / (predict.table[1, 1] + predict.table[1, 2])
treatAccuracy <- predict.table[2, 2] / (predict.table[2, 1] + predict.table[2, 2])
paste0("Con accuracy: ", sprintf("%.3f", conAccuracy))
paste0("Treat accuracy: ", sprintf("%.3f", treatAccuracy))


###=======Confusion Matrix

draw_confusion_matrix <- function(cm) {
  layout(matrix(c(1, 1, 2)))
  par(mar = c(2, 2, 2, 2))
  plot(c(100, 345), c(280, 450),
    type = "n", xlab = "", ylab = "", xaxt = "n", yaxt = "n"
  )
  title("CONFUSION MATRIX", cex.main = 1.8, col.main = "#530004")

  # create the matrix
  rect(130, 430, 230, 360,
    col = "#ea5e51",
    border = T, lwd = 2
  )
  text(178, 438, "Stable CAD", cex = 1.4, col = "grey30")

  rect(240, 430, 340, 360,
    col = "#35838d",
    border = T, lwd = 2
  )
  text(288, 438, "Acute MI", cex = 1.4, col = "grey30")

  text(100, 360, "Actual", cex = 1.5, srt = 90, font = 2, col = "#530004")
  text(235, 450, "Predicted", cex = 1.5, font = 2, col = "#530004")

  rect(130, 280, 230, 350,
    col = "#35838d",
    border = T, lwd = 2
  ) 
  rect(240, 280, 340, 350,
    col = "#ea5e51",
    border = T, lwd = 2
  ) 
  text(118, 395, "Stable CAD", cex = 1.4, srt = 90, col = "grey30")
  text(118, 315, "Acute MI", cex = 1.4, srt = 90, col = "grey30")

  # add in the cm results
  res <- as.numeric(cm$table)
  text(178, 395, res[1], cex = 1.9, font = 2, col = "white")
  text(178, 315, res[2], cex = 1.9, font = 2, col = "white")
  text(290, 395, res[3], cex = 1.9, font = 2, col = "white")
  text(290, 315, res[4], cex = 1.9, font = 2, col = "white")

  # add in the specifics
  plot(c(100, 0), c(100, 0),
    type = "n", xlab = "", ylab = "",
    main = "DETAILS", xaxt = "n", yaxt = "n", cex.main = 1.8, col.main = "#530004"
  )

  text(10, 85, names(cm$byClass[1]), cex = 1.2, font = 2, col = "#41140C")
  text(10, 65, round(as.numeric(cm$byClass[1]), 3), cex = 1.2)

  text(33, 85, names(cm$byClass[2]), cex = 1.2, font = 2, col = "#41140C")
  text(33, 65, round(as.numeric(cm$byClass[2]), 3), cex = 1.2)

  text(56, 85, names(cm$byClass[5]), cex = 1.2, font = 2, col = "#41140C")
  text(56, 65, round(as.numeric(cm$byClass[5]), 3), cex = 1.2)

  text(75, 85, names(cm$byClass[6]), cex = 1.2, font = 2, col = "#41140C")
  text(75, 65, round(as.numeric(cm$byClass[6]), 3), cex = 1.2)

  text(93, 85, names(cm$byClass[7]), cex = 1.2, font = 2, col = "#41140C")
  text(93, 65, round(as.numeric(cm$byClass[7]), 3), cex = 1.2)

  # add in the accuracy information
  text(14, 35, names(cm$byClass[3]), cex = 1.2, font = 2, col = "#41140C")
  text(14, 15, round(as.numeric(cm$byClass[6]), 3), cex = 1.2)

  text(47, 35, names(cm$byClass[4]), cex = 1.2, font = 2, col = "#41140C")
  text(47, 15, round(as.numeric(cm$byClass[7]), 3), cex = 1.2)

  text(74, 35, names(cm$overall[1]), cex = 1.2, font = 2, col = "#41140C")
  text(74, 15, round(as.numeric(cm$overall[1]), 3), cex = 1.2)

  text(94, 35, names(cm$overall[2]), cex = 1.2, font = 2, col = "#41140C")
  text(94, 15, round(as.numeric(cm$overall[2]), 3), cex = 1.2)
}

draw_confusion_matrix(confusionMatrix)


###=======ROC
y <- group
y <- ifelse(y == "CCS", 0, 1)

roc <- roc(y, as.numeric(net.predict[, 2]),
  auc = T,
  ci = T,
  percent = F
)

par(pty = "s")
plot(roc,
  main = "Test Cohort",
  print.auc = T,
  print.auc.x = 0.8,
  print.auc.y = 0.5,
  ci = TRUE,
  print.thres = TRUE,
  print.thres.cex = 1,
  legacy.axes = TRUE,
  xlab = "False Positive Percentage",
  ylab = "True Postive Percentage",
  axes = T,
  mar = c(4, 4, 2, 2) + .1,
  mgp = c(1.5, .2, 5),
  col = "#1b4965",
  lwd = 4,
  type = "l",
  MX = T,
  grid = T,
  grid.lty = 3,
  grid.lwd = 1.8,
  grid.col = "#DDDDDD",
  auc.polygon = TRUE,
  auc.polygon.col = "#91D1C266",
  identity.lty = 2,
  auc.polygon.lty = par("lty"),
  auc.polygon.density = NULL,
  auc.polygon.angle = 90,
  auc.polygon.border = T,
)
plot(ci(roc, of = "thresholds", thresholds = "best"), lwd = 2)



