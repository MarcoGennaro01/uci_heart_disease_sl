library(ucimlrepo)
library(skimr)
library(tidyr)
library(corrplot)
library(dplyr)
library(ISLR2)
library(boot) 
library(pROC)
library(MASS)
library(tree)
library(randomForest)

################################################################################
#DATA CLEANING
################################################################################

raw_data <- read.csv("./data/data.csv")
colnames(raw_data) <- c("age", "sex", "cp", "trestbps", "chol", 
                        "fbs", "restecg", "thalach", "exang", 
                        "oldpeak", "slope", "ca", "thal", "num")

data <- raw_data

data$sex <- as.factor(data$sex)
data <- data %>% rename(is_male = sex)

data$cp <- as.factor(data$cp) #chest pain type (1-4)

data$fbs <- as.factor(data$fbs) #fasting blood sugar higher than 120
data <- data %>% rename(is_diabetic = fbs)

data$restecg <- as.factor(data$restecg)

data$exang <- as.factor(data$exang)

data$slope <- as.factor(data$slope)
data$num <- as.factor(if_else(data$num==0,"0","1"))

data$thal <- factor(data$thal, levels = c("3.0", "6.0", "7.0"))
data$ca <- factor(data$ca, levels = c("0.0", "1.0", "2.0", "3.0"))

data <- drop_na(data)

################################################################################
#EDA
################################################################################

skim(data)

corrplot(cor(raw_data[, sapply(raw_data, is.numeric)]), method="number", type="lower") 

par(mfrow = c(2, 4))
for (col in names(data)[sapply(data, is.numeric)]){
  hist(data[[col]], main = col, col = "steelblue", border = "white", xlab = "")
}
par(mfrow = c(1,1))


################################################################################
#LOGISTIC REGRESSION MODEL 
################################################################################

#SUBSET SELECTION
logmodel <- glm(num~., data = data, family = "binomial") 

step_model <- step(logmodel, direction = "backward")

AIC(step_model)
AIC(logmodel)

#Weighted cost function
cost_function <- function(r, pi){
  w1 <- 4
  w0 <- 1
  cut <- 1/(1+w1/w0)
  c1 <- (r==1)&(pi < cut)
  c0 <- (r==0)&(pi >= cut)
  tc <- mean(w1*c1 + w0*c0)  
  return(tc)
}


#5-fold Cross Validation LOG MODEL
set.seed(1234)
t <- 20
kcv_err <- numeric(t)
for (i in 1:t){
  kcv_err[i] <- cv.glm(data,step_model,K=5, cost = cost_function)$delta[1]
}
(avg_kcv_err <- mean(kcv_err))
(sd(kcv_err))


#In-sample AUC LOG MODEL

predicted <- predict(step_model, data, type = "response") 
(roc_auc <- roc(data$num, predicted))


#LOOCV AUC LOG MODEL

loocv_auc_vec <- numeric(nrow(data))
for (i in 1:nrow(data)){
  tr <- data[-i,]
  te <- data[i,]
  
  model <- glm(num ~ is_male + cp + trestbps + exang + oldpeak + slope + 
                 ca + thal + thalach, data = tr, family = "binomial")
  loocv_auc_vec[i] <- predict(model,newdata = te, type = "response")
}
roc_loocv <- roc(data$num, loocv_auc_vec)
(auc(roc_loocv))


#Comparing LOOCV AUC and In-Sample AUC LOG MODEL

plot(roc_auc,  
     legacy.axes = TRUE,
     col = "steelblue",      
     lwd = 3,           
     main = "Heart Disease Model ROC Curve")

lines(roc_loocv,
      col = "darkred",
      lwd = 3)

legend("bottomright", 
       legend = c(paste("In-sample AUC =", round(auc(roc_auc), 3)),
                  paste("LOOCV AUC =", round(auc(roc_loocv), 3))),
       col = c("steelblue", "darkred"),
       lwd = 3)

#Confusion Matrix
log_weighted_classes <- as.factor(ifelse(predicted > 0.2, "1", "0"))

round(prop.table(
  table(Predicted = log_weighted_classes, Observed = data$num),
  margin = 2) * 100,
  1)

################################################################################
#LDA MODEL 
################################################################################

(lda_model <- lda(num~.,data))
lda_predicted <- predict(lda_model)

lda_scores <- lda_predicted$x[,1]
ldahist(lda_scores, g = data$num)

lda_classes <- lda_predicted$class 

round(prop.table(
  table(Predicted = lda_predicted$class, Observed = data$num),
  margin = 2) * 100,
  1)


#LOOCV AUC LDA
loocv_lda_vec <- numeric(nrow(data))

for (i in 1:nrow(data)){
  tr <- data[-i, ]
  te <- data[i, ]
  
  model <- lda(num ~ ., data = tr)
  loocv_lda_vec[i] <- predict(model, newdata = te)$posterior[, 2]
}

roc_lda_loocv <- roc(data$num, loocv_lda_vec)
(auc(roc_lda_loocv))


#5-fold Cross Validation LDA

set.seed(1234)
t <- 20
kcv_lda_err <- numeric(t)

for (i in 1:t){
  k <- 5
  folds <- sample(rep(1:k, length.out = nrow(data)))
  fold_costs <- numeric(k)
  
  for (j in 1:k){
    tr <- data[folds != j, ]
    te <- data[folds == j, ]
    
    model <- lda(num ~ ., data = tr)
    preds <- predict(model, newdata = te)$posterior[, 2]
    fold_costs[j] <- cost_function(as.numeric(as.character(te$num)), preds)
  }
  kcv_lda_err[i] <- mean(fold_costs)
}

(avg_kcv_lda_err <- mean(kcv_lda_err))
(sd(kcv_lda_err))


#Confusion Matrix with the same threshold 0.2

lda_weighted_classes <- as.factor(ifelse(lda_predicted$posterior[, 2] > 0.2, 
                                         "1",
                                         "0"))

conf_matrix_lda <- table(Predicted = lda_weighted_classes, Actual = data$num)
round(prop.table(conf_matrix_lda, margin = 2) * 100, 1)


################################################################################
#TREE BASED METHODS
################################################################################

tree_model <- tree(num ~ ., data = data)
cv_tree_model <- cv.tree(tree_model, FUN = prune.misclass)

plot(tree_model, 
     lwd = 1.5,           
     type = "uniform",    
     col = "black")   

text(tree_model, 
     cex = 1,          
     col = "black",       
     pretty = 0,          
     font = 2)            

par(mfrow = c(1, 2))
plot(cv_tree_model$size,cv_tree_model$dev,type="b",lwd=2,col="blue",
     xlab="Number of terminal nodes", ylab="Deviance" )
plot(cv_tree_model$k,cv_tree_model$dev,type="b",lwd=2,col="blue",
     xlab="Cost-complexity", ylab="Deviance" )
par(mfrow = c(1, 1))

#Pruning
size_pruned_tree <- cv_tree_model$size[which.min(cv_tree_model$dev)]

pruned_tree <- prune.misclass(tree_model,best=size_pruned_tree)

plot(pruned_tree, lwd = 1.5, type = "uniform", col = "black")
text(pruned_tree, cex = 1, col = "black", pretty = 0, font = 2)
title(main = "Pruned Decision Tree - Heart Disease")

#LOOCV AUC PRUNED TREE
loocv_tree_vec <- numeric(nrow(data))

for (i in 1:nrow(data)){
  tr <- data[-i, ]
  te <- data[i, ]
  
  model <- tree(num ~ ., data = tr)
  pruned <- prune.misclass(model, best = size_pruned_tree)
  loocv_tree_vec[i] <- predict(pruned, newdata = te, type = "vector")[, 2]
}

roc_tree_loocv <- roc(data$num, loocv_tree_vec)
(auc(roc_tree_loocv))


#5-fold Cross Validation PRUNED TREE
set.seed(1234)
t <- 20
kcv_tree_err <- numeric(t)

for (i in 1:t){
  k <- 5
  folds <- sample(rep(1:k, length.out = nrow(data)))
  fold_costs <- numeric(k)
  
  for (j in 1:k){
    tr <- data[folds != j, ]
    te <- data[folds == j, ]
    
    model <- tree(num ~ ., data = tr)
    pruned <- prune.misclass(model, best = size_pruned_tree)
    preds <- predict(pruned, newdata = te, type = "vector")[, 2]
    fold_costs[j] <- cost_function(as.numeric(as.character(te$num)), preds)
  }
  kcv_tree_err[i] <- mean(fold_costs)
}

(avg_kcv_tree_err <- mean(kcv_tree_err))
(sd(kcv_tree_err))

################################################################################
#RANDOM FOREST
################################################################################

set.seed(1234)
(rf_model <- randomForest(num ~ ., 
                          mtry = ncol(data)-1,
                          data = data,
                          importance = TRUE))

#LOOCV AUC RANDOM FOREST
loocv_rf_vec <- numeric(nrow(data))

for (i in 1:nrow(data)){
  tr <- data[-i, ]
  te <- data[i, ]
  
  model <- randomForest(num ~ ., data = tr, mtry = ncol(data)-1)
  loocv_rf_vec[i] <- predict(model, newdata = te, type = "prob")[, 2]
}

roc_rf_loocv <- roc(data$num, loocv_rf_vec)
(auc(roc_rf_loocv))

#5-fold Cross Validation RANDOM FOREST
set.seed(1234)
t <- 20
kcv_rf_err <- numeric(t)

for (i in 1:t){
  k <- 5
  folds <- sample(rep(1:k, length.out = nrow(data)))
  fold_costs <- numeric(k)
  
  for (j in 1:k){
    tr <- data[folds != j, ]
    te <- data[folds == j, ]
    
    model <- randomForest(num ~ ., data = tr, mtry = ncol(data)-1)
    preds <- predict(model, newdata = te, type = "prob")[, 2]
    fold_costs[j] <- cost_function(as.numeric(as.character(te$num)), preds)
  }
  kcv_rf_err[i] <- mean(fold_costs)
}

(avg_kcv_rf_err <- mean(kcv_rf_err))
(sd(kcv_rf_err))

varImpPlot(rf_model, col="steelblue", pch=25, main="Importance Plot")

################################################################################
#FINAL MODEL COMPARISON
################################################################################

comparison_table <- data.frame(
  Model = c("Step Model", "LDA", "Pruned Tree", "Random Forest"),
  LOOCV_AUC = c(round(auc(roc_loocv), 3),
                round(auc(roc_lda_loocv), 3),
                round(auc(roc_tree_loocv), 3),
                round(auc(roc_rf_loocv), 3)),
  CV_Error_Mean = round(c(avg_kcv_err, avg_kcv_lda_err, 
                          avg_kcv_tree_err, avg_kcv_rf_err), 3),
  CV_Error_SD = round(c(sd(kcv_err), sd(kcv_lda_err), 
                        sd(kcv_tree_err), sd(kcv_rf_err)), 3)
)

print(comparison_table)
