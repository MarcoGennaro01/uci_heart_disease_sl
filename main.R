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

data <- data %>% mutate(across(where(is.numeric), scale))

################################################################################
#EDA
################################################################################

skim(data)

corrplot(cor(drop_na(raw_data)), method="number", type="lower") 

################################################################################
#REGRESSION MODEL 
################################################################################

#Logistic regression
logmodel <- glm(num~., data = data, family = "binomial") 

#Defining the weighted cost function
cost_function <- function(r, pi){
  w1 <- 4
  w0 <- 1
  cut <- 1/(1+w1/w0)
  c1 <- (r==1)&(pi < cut)
  c0 <- (r==0)&(pi >= cut)
  tc <- mean(w1*c1 + w0*c0)  
  return(tc)
}

#5-fold Cross Validation
set.seed(1234)
t <- 20
kcv_err <- numeric(t)
for (i in 1:t){
  kcv_err[i] <- cv.glm(data,logmodel,K=5, cost = cost_function)$delta[1]
}
(avg_kcv_err <- mean(kcv_err))
(sd(kcv_err))

#AUC 
predicted <- predict(logmodel, data, type = "response") 
roc_object <- roc(data$num, predicted)
cat("Area under the curve:", auc(roc_object))
plot(roc_object, 
     legacy.axes = TRUE,
     col = "blue",      
     lwd = 3,           
     main = "Heart Disease Model ROC Curve")


#LDA
lda_model <- lda(num~.,data)
lda_model

