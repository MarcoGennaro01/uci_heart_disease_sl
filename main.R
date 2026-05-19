library(ucimlrepo)
library(skimr)
library(tidyr)
library(corrplot)
library(dplyr)
library(ISLR2)
library(boot) 
library(pROC)
library(MASS)

################################################################################
#DATA CLEANING
################################################################################

raw_data <- fetch_ucirepo(id = 45)$data$original 
data <- drop_na(raw_data)

data$sex <- as.factor(data$sex)
data <- data %>% rename(is_male = sex)

data$cp <- as.factor(data$cp) #chestpain type (1-4)

data$fbs <- as.factor(data$fbs) #fasting blood sugar higher than 120
data <- data %>% rename(is_diabetic = fbs)

data$restecg <- as.factor(data$restecg)

data$exang <- as.factor(data$exang)

data$slope <- as.factor(data$slope)
data$num <- as.factor(if_else(data$num==0,"0","1"))


################################################################################
#EDA
################################################################################

skim(data)

corrplot(cor(drop_na(raw_data)), method="number", type="lower") 

par(mfrow = c(2, 4))
for (col in names(data)[sapply(data, is.numeric)]){
  hist(data[[col]], main = col, col = "steelblue", border = "white", xlab = "")
}
par(mfrow = c(1,1))

################################################################################
#LOGISTIC REGRESSION MODEL 
################################################################################

logmodel <- glm(num~., data = data, family = "binomial") 


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
  kcv_err[i] <- cv.glm(data,logmodel,K=5, cost = cost_function)$delta[1]
}
(avg_kcv_err <- mean(kcv_err))
(sd(kcv_err))


#In-sample AUC LOG MODEL

predicted <- predict(logmodel, data, type = "response") 
(roc_auc <- roc(data$num, predicted))
(auc(roc_auc))


#LOOCV AUC LOG MODEL

loocv_auc_vec <- numeric(nrow(data))
for (i in 1:nrow(data)){
  tr <- data[-i,]
  te <- data[i,]
  
  model <- glm(num~., data = tr, family = "binomial")
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
round(prop.table(
  table(Predicted = log_classes, Observed = data$num),
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
table(Predicted = lda_predicted$class, Observed = data$num)


#LOOCV AUC LCA
loocv_lda_vec <- numeric(nrow(data))

for (i in 1:nrow(data)){
  tr <- data[-i, ]
  te <- data[i, ]
  
  model <- lda(num ~ ., data = tr)
  loocv_lda_vec[i] <- predict(model, newdata = te)$posterior[, 2]
}

roc_lda_loocv <- roc(data$num, loocv_lda_vec)
(auc(roc_lda_loocv))


#5-fold Cross Validation LCA

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
conf_matrix_lda <- table(Predicted = lda_classes_thresh, Actual = data$num)
round(prop.table(conf_matrix_lda, margin = 2) * 100, 1)

