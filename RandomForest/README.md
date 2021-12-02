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
```
