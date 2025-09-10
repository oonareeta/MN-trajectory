# Libraries
source("mounts/research/src/Rfunctions/library.R")
library(caret)
library(randomForestSRC)
library(tidyr)
library(pROC)
library(MLmetrics)
library(xgboost)
library(doParallel)
library(PRROC)
library(cutpointr)
library(SHAPforxgboost)

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")



### AML

# Load the model
AML_model <- xgb.load("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/de_novo_AML_basic_model.json")

# Load x_train
x_train = fread("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/de_novo_AML_x_train.csv")
x_train = as.matrix(x_train)

# SHAP
shap_values <- shap.values(xgb_model = AML_model, X_train = x_train)
tmp4 = shap_values$shap_score
fwrite(tmp4, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/de_novo_AML_basic_model_feature_values.csv")


tmp3 = data.frame(names(shap_values$mean_shap_score), shap_values$mean_shap_score)
colnames(tmp3)[1] = "names";colnames(tmp3)[2] = "mean_shap_score"
fwrite(tmp3, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/de_novo_AML_basic_model_SHAP_scores.csv")

### MDS

# Load the model
MDS_model <- xgb.load("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/MDS_basic_model.json")

# Load x_train
x_train = fread("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/MDS_x_train.csv")
x_train = as.matrix(x_train)

# SHAP
shap_values <- shap.values(xgb_model = MDS_model, X_train = x_train)
tmp4 = shap_values$shap_score
fwrite(tmp4, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/MDS_basic_model_feature_values.csv")


tmp3 = data.frame(names(shap_values$mean_shap_score), shap_values$mean_shap_score)
colnames(tmp3)[1] = "names";colnames(tmp3)[2] = "mean_shap_score"
fwrite(tmp3, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/MDS_basic_model_SHAP_scores.csv")


### MF

# Load the model
MF_model <- xgb.load("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/MF_basic_model.json")

# Load x_train
x_train = fread("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/MF_x_train.csv")
x_train = as.matrix(x_train)

# SHAP
shap_values <- shap.values(xgb_model = MF_model, X_train = x_train)
tmp4 = shap_values$shap_score
fwrite(tmp4, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/MF_basic_model_feature_values.csv")


tmp3 = data.frame(names(shap_values$mean_shap_score), shap_values$mean_shap_score)
colnames(tmp3)[1] = "names";colnames(tmp3)[2] = "mean_shap_score"
fwrite(tmp3, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/MF_basic_model_SHAP_scores.csv")



### any_MN

# Load the model
any_MN_model <- xgb.load("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/any_MN_basic_model.json")

# Load x_train
x_train = fread("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/any_MN_x_train.csv")
x_train = as.matrix(x_train)

# SHAP
shap_values <- shap.values(xgb_model = any_MN_model, X_train = x_train)
tmp4 = shap_values$shap_score
fwrite(tmp4, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/any_MN_basic_model_feature_values.csv")


tmp3 = data.frame(names(shap_values$mean_shap_score), shap_values$mean_shap_score)
colnames(tmp3)[1] = "names";colnames(tmp3)[2] = "mean_shap_score"
fwrite(tmp3, "mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/any_MN_basic_model_SHAP_scores.csv")
