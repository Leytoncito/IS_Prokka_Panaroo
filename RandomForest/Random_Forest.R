#import libraries

library(randomForest)
library(ggplot2)
library(caret)
library(pROC)
library(readxl)
library(caTools)

#import metadata

df<-read_xlsx("Supplementary_Data_1.xlsx", sheet = "detailed_dataframe")
df_lineage = df[, c(2, 22:351)]
df_location = df[, c(3, 22:351)]
df_arg = df[, c(6, 22:351)]
df_lineage$Lineage<-as.factor(df_lineage$Lineage)
df_location$Location<-as.factor(df_location$Location)
df_arg$`No ARG`<-as.factor(df_arg$`No ARG`)

#divide set into testing set and training set

set.seed(5)
split = sample.split(df_location$Location, SplitRatio = 0.80)
training_set = subset(df_location, split == TRUE)
testing_set = subset(df_location, split == FALSE)

set.seed(5)
split = sample.split(df_lineage$Lineage, SplitRatio = 0.80)
training_set = subset(df_lineage, split == TRUE)
testing_set = subset(df_lineage, split == FALSE)

set.seed(5)
split = sample.split(df_arg$`No ARG`, SplitRatio = 0.80)
training_set = subset(df_arg, split == TRUE)
testing_set = subset(df_arg, split == FALSE)


# RF models

rf_model_lineage <- randomForest(Lineage ~., data = training_set,
                                ntree=501, importance=T, confusion=T, err.rate=T)


rf_model_location <- randomForest(Lineage ~., data = training_set,
                                ntree=501, importance=T, confusion=T, err.rate=T)

rf_model_arg <- randomForest(Lineage ~., data = training_set,
                                ntree=501, importance=T, confusion=T, err.rate=T)


#Kappa test 
pred_test_class1 <- predict(rf_model_lineage, testing_set,
                            type="class", norm.votes=TRUE, predict.all=FALSE,
                            proximity=FALSE, nodes=FALSE)
head(pred_test_class1)
Cm<-confusionMatrix(pred_test_class, testing_set$Lineage)
table1<-Cm$table

pred_test_class2 <- predict(rf_model_location, testing_set,
                            type="class", norm.votes=TRUE, predict.all=FALSE,
                            proximity=FALSE, nodes=FALSE)
head(pred_test_class2)
Cm<-confusionMatrix(pred_test_class, testing_set$Lineage)
table2<-Cm$table

pred_test_class3 <- predict(rf_model_arg, testing_set,
                            type="class", norm.votes=TRUE,
                            predict.all=FALSE, proximity=FALSE, nodes=FALSE)
head(pred_test_class3)
Cm<-confusionMatrix(pred_test_class, testing_set$Lineage)
table3<-Cm$table

library(fmsb)
kappa_lineage<-Kappa.test(table1)
kappa_lineage

kappa_location<-Kappa.test(table2)
kappa_location

kappa_arg<-Kappa.test(table3)
kappa_arg

# Cross-Validation

library(rfUtilities)
fold_lineage<-rf.crossValidation(rf_model_lineage, training_set, ydata = NULL, p = 0.2, n = 100,
                                seed = 1, normalize = FALSE, bootstrap = FALSE, trace = TRUE)

fold_location<-rf.crossValidation(rf_model_location, training_set, ydata = NULL, p = 0.2, n = 100,
                                seed = 1, normalize = FALSE, bootstrap = FALSE, trace = TRUE)

fold_arg<-rf.crossValidation(rf_model_arg, training_set, ydata = NULL, p = 0.2, n = 100,
                                seed = 1, normalize = FALSE, bootstrap = FALSE, trace = TRUE)

#save output
library(writexl)
write_xlsx(rf_model_arg$importance, file = "arg_gini.xlsx")
write_xlsx(rf_model_linaje$importance, file= "lineage_gini.xlsx")
write_xlsx(rf_model_location$importance, file="location_gini.xlsx")