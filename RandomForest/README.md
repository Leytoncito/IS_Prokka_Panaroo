## Random forest (Summary code)

First, We need import necessary libraries.

```
library(randomForest)
require(pROC)
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

We divide the metadata in testing set and training set. We also established seeds for reproduction

```
set.seed(5) # predice todos los linajes
split = sample.split(df_lineage$Lineage, SplitRatio = 0.80)
training_set = subset(df_lineage, split == TRUE)
testing_set = subset(df_lineage, split == FALSE)

```

Then we create the RF model

```
rf_model_linaje <- randomForest(Lineage ~., data = training_set,
                                ntree=501, importance=T, confusion=T, err.rate=T)
```

We validate the model with the test set, perform a kappa test, and a Cross-Validation.
                                
```
pred_test_class <- predict(rf_model_linaje, testing_set, type="class", norm.votes=TRUE, predict.all=FALSE, proximity=FALSE, nodes=FALSE)
head(pred_test_class)
Cm<-confusionMatrix(pred_test_class, testing_set$Lineage)
table<-Cm$table

##Kappa test
kappa_lineage<-Kappa.test(table)
kappa_lineage

#Cross-Validation
fold_linaje<-rf.crossValidation(rf_model_linaje, training_set, ydata = NULL, p = 0.1, n = 100,
                   seed = 5, normalize = FALSE, bootstrap = FALSE, trace = FALSE)

```
Finally, we save the files with gini index.
```
write.xlsx(rf_model_linaje$importance, file= "lineage_gini.xlsx")
```

