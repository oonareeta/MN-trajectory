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

disease = "MDS"

# Load the model
MDS_model <- xgb.load(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_git/results/final_model/", disease, "_final_model.json"))

# Load x_train
x_train = read.csv(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_git/results/final_model/SHAP/", disease, "_x_train.csv"))
x_train = as.matrix(x_train)

# SHAP
shap_values <- shap.values(xgb_model = MDS_model, X_train = x_train) 

tmp4 = shap_values$shap_score
fwrite(tmp4, paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_git/results/final_model/SHAP/", disease, "_final_model_feature_values.csv"))

tmp3 = data.frame(names(shap_values$mean_shap_score), shap_values$mean_shap_score)
colnames(tmp3)[1] = "names";colnames(tmp3)[2] = "mean_shap_score"
fwrite(tmp3, paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_git/results/final_model/SHAP/", disease, "_final_model_SHAP_scores.csv"))




