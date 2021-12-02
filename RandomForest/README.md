## Random forest

First, We need import necessary libraries.

```
library(randomForest)
require(pROC)
require(raster)
require(rgdal)
require(tmap)
require(ggplot2)
require(caret)
library(readxl)
library(caTools)
library(fmsb)
library(rfUtilities)
```

Then, we import the dataframe and we created subset

```
df<-read_xlsx("Supplementary_Data_1.xlsx", sheet = "detailed_dataframe")
df_lineage = df[, c(2, 22:351)]
df_location = df[, c(3, 22:351)]
df_arg = df[, c(6, 22:351)]
df_lineage$Lineage<-as.factor(df_lineage$Lineage)
df_location$Location<-as.factor(df_location$Location)
df_arg$`No ARG`<-as.factor(df_arg$`No ARG`)
```
```
set.seed(5) # predice todos los linajes
set.seed(42) # significacion de location es 0,06.
set.seed(1234)
set.seed(34) #Decente para location a 0.8
set.seed(36)# p-value 0.02 para location 0.8
split = sample.split(df_lineage$Lineage, SplitRatio = 0.80)
training_set = subset(df_lineage, split == TRUE)
testing_set = subset(df_lineage, split == FALSE)

rf_model_linaje <- randomForest(Lineage ~., data = training_set,
                                ntree=501, importance=T, confusion=T, err.rate=T)
                                
pred_test_class <- predict(rf_model_linaje, testing_set, type="class", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
head(pred_test_class)
Cm<-confusionMatrix(pred_test_class, testing_set$Lineage)
table<-Cm$table

library(fmsb)
kappa_location<-Kappa.test(table)
kappa_location
kappa_lineage<-Kappa.test(table)
kappa_lineage
kappa_arg<-Kappa.test(table)
kappa_arg

write.xlsx(rf_model_arg$importance, file = "arg_gini.xlsx")
write.xlsx(rf_model_linaje$importance, file= "lineage_gini.xlsx")
write.xlsx(rf_model_location$importance, file="location_gini.xlsx")

#Cross-Validation
library(rfUtilities)
fold_linaje<-rf.crossValidation(rf_model_linaje, training_set, ydata = NULL, p = 0.5, n = 100,
                   seed = NULL, normalize = FALSE, bootstrap = FALSE, trace = FALSE)

significance_linaje<-rf.significance(rf_model_linaje, training_set[-1], nperm = 10, p= 0.05)
significance_linaje
```
