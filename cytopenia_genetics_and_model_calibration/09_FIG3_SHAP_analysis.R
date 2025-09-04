# Predict AML, MDS and MF using XGBoost

# Libraries
source("mounts/research/src/Rfunctions/library.R")
library(tidyr)
library(xgboost)
library(SHAPforxgboost)

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")


# Load data
for (i in c("any_MN", "MDS", "MF", "de novo AML")) {
  print(paste0("Processing ", i))
  
  if (i == "de novo AML") {
    j = "de_novo_AML"
  } else {
    j = i
  }
  
  # Load the model
  model1 <- xgb.load(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/", j, "_basic_model.json"))
  
  # Load x_train
  x_train = fread(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/results/basic_model/SHAP/", j, "_x_train.csv"), nrows=500000)
  x_train = as.matrix(x_train)
  
  # SHAP
  shap_values <- shap.values(xgb_model = model1, X_train = x_train)
  
  # The ranked features by mean |SHAP|
  # To prepare the long-format data:
  shap_long <- shap.prep(shap_contrib = shap_values$shap_score, X_train = x_train) %>%# as.matrix(train_data[,5:(ncol(train_data)-1)]))
    dplyr::filter(variable %in% names(shap_values$mean_shap_score[1:10]))
  
  dir.create(paste0(results, "/", i, "/SHAP_final"))
  
  # Correct names
  names1 = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/variable_names.csv")
  shap_long = shap_long %>%
    dplyr::left_join(names1) %>%
    dplyr::rename(variable1 = output)
  
  for ( k1 in 1:10) {
    
    print(names(shap_values$mean_shap_score)[k1])
    
    # Select variable
    var1 = shap_long %>%
      dplyr::filter(variable == names(shap_values$mean_shap_score)[k1]) %>%
      distinct(variable1)
    
    # Plot
    g = shap.plot.dependence(data_long = shap_long %>%
                               dplyr::filter(variable == names(shap_values$mean_shap_score)[k1]) %>%
                               slice_sample(n=10000), names(shap_values$mean_shap_score)[k1]) +
      theme_bw() +
      xlab(var1$variable1) +
      ylab("SHAP value") +
      theme(axis.text.x = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title = element_text(size=14, colour = "black"),
            axis.line = element_line(colour = "black"),
            panel.grid.major = element_blank(),
            panel.grid.minor = element_blank(),
            panel.border = element_blank(),
            panel.background = element_blank()); g
    
    ggsave(plot = g, filename = paste0(results, "/", i, "/SHAP_final/SHAP_contribution_", names(shap_values$mean_shap_score)[k1], ".png"), width = 5, height = 5, units = "in", dpi = 300)
  }
  
  gc()
  
}
