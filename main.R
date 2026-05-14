library(ucimlrepo)
library(skimr)
library(tidyr)
library(corrplot)
library(dplyr)
library(ISLR2)

################################################################################
#DATA CLEANING
################################################################################

raw_data <- fetch_ucirepo(id = 45)$data$original 
data <- drop_na(raw_data)

data$sex=as.factor(data$sex)
data <- data %>% rename(is_male = sex)

data$cp=as.factor(data$cp) #chestpain type (1-4)

data$fbs=as.factor(data$fbs) #fasting blood sugar higher than 120
data <- data %>% rename(is_diabetic = fbs)

data$restecg=as.factor(data$restecg)

data$exang=as.factor(data$exang)

data$slope=as.factor(data$slope)

data$num=as.factor(if_else(data$num==0,"0","1"))

data <- data %>% mutate(across(where(is.numeric), scale))

################################################################################
#EDA
################################################################################

skim(data)

corrplot(cor(raw_data), method="number", type="lower") 

################################################################################
#REGRESSION
################################################################################

logmodel <- glm(num~., data = data, family = "binomial") 

predicted <- predict.glm(logmodel,type="response")

predicted=(if_else(predicted>0.5,"1","0"))
table(predicted, data$num)
table(glm.pred2,Default$default) # confusion matrix









