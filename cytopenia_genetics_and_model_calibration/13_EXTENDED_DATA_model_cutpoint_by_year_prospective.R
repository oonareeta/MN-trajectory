# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Libraries
source("mounts/research/src/Rfunctions/library.R")
library(caret)
library(pROC)
library(MLmetrics)
library(PRROC)


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Load test data
disease = c("any_MN", "de_novo_AML", "MDS", "MF")
metrics1 = data.frame()

# Loop over diseases
for (i in disease) {
  
  for (j in c(1,3,5)) {
  
  # Load test data
  df = read.csv(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/trajectory_model/XGBoost_COX/youden/", i, "_test_data_with_predictions.csv"))
  
  # Replace min time_to_dg as 0
  df_tmp = df %>%
    mutate(time_to_dg = abs(time_to_dg)) %>%
    group_by(patient_id) %>%
    summarise(fu_time = max(time_to_dg, na.rm=TRUE)) %>%
    ungroup()
  
  # Select only disease patients
  df_disease = df %>%
    mutate(time_to_dg = abs(time_to_dg)) %>%
    dplyr::left_join(df_tmp) %>%
    dplyr::mutate(fu_time = 1+fu_time-time_to_dg) %>%
    dplyr::filter(fu_time <= j*365.24) %>%
    dplyr::rename(disease = disease_status) %>%
    dplyr::mutate(risk_score = log(risk_score, 10),
                  time_to_dg = abs(time_to_dg),
                  time_to_dg_year = -time_to_dg / 365.24)
  
  # Number of samples per year
  tt = df_disease %>%
    group_by(patient_id) %>%
    mutate(n = n()) %>%
    arrange(desc(fu_time)) %>%
    slice(1) %>%
    ungroup() %>%
    mutate(n1 = 365.24*n/fu_time)
  
  # Calculate the AUC
  auc = pROC::roc(df_disease$disease, df_disease$risk_score, plot=TRUE);auc
  
  # Confusion matrix
  conf_matrix <- confusionMatrix(as.factor(df_disease$predicted_label), as.factor(df_disease$disease), positive = "1")
  conf_matrix1 = as.data.frame(conf_matrix$table) %>%
    group_by(Reference) %>%
    dplyr::mutate(Prop = 100*Freq / sum(Freq)) %>%
    ungroup() %>%
    dplyr::mutate(Prediction = factor(Prediction, levels = c(0, 1)),
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
         filename = paste0(results, "/", gsub("de_novo_AML", "de novo AML", i), "/Prediction_with_max_", j, "_year_of_followup.png"),
         height = 3.5, width = 3.5, bg = "white", units = "in", dpi = 300)
  
  
  # Combine metrics
  TP = conf_matrix$table[2, 2]
  FN = conf_matrix$table[1, 2]
  TN = conf_matrix$table[1, 1]
  FP = conf_matrix$table[2, 1]
  
  # Print the results
  metrics = data.frame(
    disease = i,
    year = j,
    n_labtests_per_year = median(tt$n1),
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
writexl::write_xlsx(metrics1, paste0(results, "/metrics_prediction_with_max_1_3_5_years_of_followup.xlsx"))
