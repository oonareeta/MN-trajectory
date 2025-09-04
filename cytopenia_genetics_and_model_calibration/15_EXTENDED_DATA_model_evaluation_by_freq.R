# How sampling freq is associated with performance

# Libraries
source("mounts/research/src/Rfunctions/library.R")
library(caret)
library(pROC)

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Load test data
disease = c("any_MN", "de_novo_AML", "MDS", "MF")
metrics1 = data.frame()

# Loop over diseases
for (i in disease) {
  
  # Load test data
  df = read.csv(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/trajectory_model/XGBoost_COX/youden/", i, "_test_data_with_predictions.csv"))
  
  # Select only disease patients
  df_disease = df %>%
    # dplyr::select(-risk_score) %>%
    dplyr::rename(disease = disease_status) %>%
    dplyr::mutate(risk_score = log(risk_score, 10),
                  time_to_dg_year = -time_to_dg / 365.24)
  
  # Folder
  # if (i == "primary_MF") {
  #   j = "MF"
  if (i == "de_novo_AML") {
    j = "de novo AML"
  } else {
    j = i
  }
  
  # Divide samples into 5 groups
  n = 5
  df_disease$n_prev_3 = ntile(df_disease$n_prev, n)
  for (k in 1:n) {
    df_disease1 = df_disease[df_disease$n_prev_3==k,]
    k1 = min(df_disease1$n_prev)
    
    # Calculate the AUC
    auc = pROC::roc(df_disease1$disease, df_disease1$risk_score, plot=TRUE);print(auc)
    
    # Confusion matrix
    conf_matrix <- confusionMatrix(as.factor(df_disease1$predicted_label), as.factor(df_disease1$disease), positive = "1")
    conf_matrix1 = as.data.frame(conf_matrix$table) %>%
      group_by(Reference) %>%
      dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
      ungroup() %>%
      dplyr::mutate(
        Prediction = factor(Prediction, levels = c(0, 1)),
        Reference = factor(Reference, levels = c(0, 1)))
    
    # Color
    if (i == "any_MN") {
      high = "#2D0E3D"
    } else if (i == "MDS") {
      high = "#348ABD"
    } else if (i == "MF") {
      high = "#2B6E2A"
    } else if (i == "de_novo_AML") {
      high = "#BF9F45"
    }
    
    # Identify the top-right (FP) and bottom-left (FN) cells
    conf_matrix1 <- conf_matrix1 %>%
      mutate(
        FontColor = case_when(
          (Prop>50 | (high == "#2D0E3D" & Prop>40)) ~ "white",  # Top-right (FP)
          TRUE ~ "black"
        )
      )
    
    # Plot
    # Create a nice plot using ggplot2
    g = ggplot(data = conf_matrix1, aes(y = reorder(Reference, -as.integer(Reference)),
                                        x = reorder(Prediction, as.integer(Prediction)), fill = Prop)) +
      geom_tile(color="black", linewidth = 0.6) +
      geom_text(aes(label=paste0(round(Freq, 0), "\n", ifelse(str_detect(round(Prop, 1), "\\."), round(Prop, 1), paste0(round(Prop, 1), ".0")), "%"),
                    color = FontColor), size=4.5) +
      # scale_fill_distiller(palette = "Reds", direction = 1) +
      scale_fill_gradient(low="white", high=high) +
      scale_color_identity() +
      labs(x = "Predicted label", y = "True label", fill = "Proportion (%)") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_text(colour = "black", size = 12),
            axis.text.y = element_text(colour = "black", size = 11),
            axis.text.x = element_text(colour = "black", size = 11),
            legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); g
    ggsave(plot = g,
           filename = paste0(results, "/", gsub("de_novo_AML", "de novo AML", i), "/Prediction_with_min_", k1, "_previous_labs.png"),
           height = 3.5, width = 3.5, bg = "white", units = "in", dpi = 300)
    
    
    # Combine metrics
    TP = conf_matrix$table[2, 2]
    FN = conf_matrix$table[1, 2]
    TN = conf_matrix$table[1, 1]
    FP = conf_matrix$table[2, 1]
    
    # Print the results
    metrics = data.frame(
      disease = i,
      n_prev_samples = k1,
      TP = TP,
      FP = FP,
      TN = TN,
      FN = FN,
      auc = round(auc$auc, 3),
      precision = round(TP / (TP + FP), 3),
      recall = round(TP / (TP + FN), 3)
    ) %>%
      dplyr::mutate(f1 = round((2*precision*recall)/(precision+recall), 3))
    
    # Rbind
    metrics1 = rbind(metrics1, metrics)
    
  }
}

# Export
writexl::write_xlsx(metrics1, paste0(results, "/metrics_prediction_n_samples.xlsx"))


############### The same but focusing on the last 2 years ###############


metrics1 = data.frame()

# Loop over diseases
for (i in disease) {
  
  # Load test data
  df = read.csv(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/trajectory_model/XGBoost_COX/youden/", i, "_test_data_with_predictions.csv"))
  
  # Select only disease patients
  df_disease = df %>%
    # dplyr::select(-risk_score) %>%
    dplyr::rename(disease = disease_status) %>%
    dplyr::mutate(risk_score = log(risk_score, 10),
                  time_to_dg = abs(time_to_dg),
                  time_to_dg_year = -time_to_dg / 365.24) %>%
    dplyr::filter(time_to_dg <= (2*365.24))
  
  # Folder
  if (i == "de_novo_AML") {
    j = "de novo AML"
  } else {
    j = i
  }
  
  n = 5
  df_disease$n_prev_3 = ntile(df_disease$n_prev, n)
  for (k in 1:5) {
    df_disease1 = df_disease[df_disease$n_prev_3==k,]
    k1 = min(df_disease1$n_prev)
    
    # Calculate the AUC
    auc = pROC::roc(df_disease1$disease, df_disease1$risk_score, plot=TRUE);print(auc)
    
    # Confusion matrix
    conf_matrix <- confusionMatrix(as.factor(df_disease1$predicted_label), as.factor(df_disease1$disease), positive = "1")
    conf_matrix1 = as.data.frame(conf_matrix$table) %>%
      group_by(Reference) %>%
      dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
      ungroup() %>%
      dplyr::mutate(
        Prediction = factor(Prediction, levels = c(0, 1)),
        Reference = factor(Reference, levels = c(0, 1)))
    
    # Color
    if (i == "any_MN") {
      high = "#2D0E3D"
    } else if (i == "MDS") {
      high = "#348ABD"
    } else if (i == "MF") {
      high = "#2B6E2A"
    } else if (i == "de_novo_AML") {
      high = "#BF9F45"
    }
    
    # Identify the top-right (FP) and bottom-left (FN) cells
    conf_matrix1 <- conf_matrix1 %>%
      mutate(
        FontColor = case_when(
          (Prop>50 | (high == "#2D0E3D" & Prop>40)) ~ "white",  # Top-right (FP)
          TRUE ~ "black"
        )
      )
    
    # Plot
    # Create a nice plot using ggplot2
    g = ggplot(data = conf_matrix1, aes(y = reorder(Reference, -as.integer(Reference)),
                                        x = reorder(Prediction, as.integer(Prediction)), fill = Prop)) +
      geom_tile(color="black", linewidth = 0.6) +
      geom_text(aes(label=paste0(round(Freq, 0), "\n", ifelse(str_detect(round(Prop, 1), "\\."), round(Prop, 1), paste0(round(Prop, 1), ".0")), "%"),
                    color = FontColor), size=4.5) +
      scale_fill_gradient(low="white", high=high) +
      scale_color_identity() +
      labs(x = "Predicted label", y = "True label", fill = "Proportion (%)") +
      theme_minimal() +
      theme(legend.position = "none",
            axis.title = element_text(colour = "black", size = 12),
            axis.text.y = element_text(colour = "black", size = 11),
            axis.text.x = element_text(colour = "black", size = 11),
            legend.title = element_text(colour = "black", size = 12, vjust = 0.8)); g
    ggsave(plot = g,
           filename = paste0(results, "/", gsub("de_novo_AML", "de novo AML", i), "/Prediction_last_2_years_with_min_", k1, "_previous_labs.png"),
           height = 3.5, width = 3.5, bg = "white", units = "in", dpi = 300)
    
    
    # Combine metrics
    TP = conf_matrix$table[2, 2]
    FN = conf_matrix$table[1, 2]
    TN = conf_matrix$table[1, 1]
    FP = conf_matrix$table[2, 1]
    
    # Print the results
    metrics = data.frame(
      disease = i,
      n_prev_samples = k1,
      TP = TP,
      FP = FP,
      TN = TN,
      FN = FN,
      auc = round(auc$auc, 3),
      precision = round(TP / (TP + FP), 3),
      recall = round(TP / (TP + FN), 3)
    ) %>%
      dplyr::mutate(f1 = round((2*precision*recall)/(precision+recall), 3))
    
    # Rbind
    metrics1 = rbind(metrics1, metrics)
    
  }
}

# Export
writexl::write_xlsx(metrics1, paste0(results, "/metrics_prediction_last_2_years_n_samples_last_2y.xlsx"))

