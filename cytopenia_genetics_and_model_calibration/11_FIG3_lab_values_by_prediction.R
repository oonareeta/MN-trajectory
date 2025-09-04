# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Libraries
source("mounts/research/src/Rfunctions/library.R")


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

## Healthy
### Loop
filelist = list.files(path = "mounts/research/husdatalake/disease/data/Preleukemia/multilab3", pattern = "file[[:digit:]]*\\.parquet", full.names = TRUE)
healthy = data.frame()
for (i in 1:10) {
  
  print(paste0("Processing ", i))
  
  # Read demo
  i=1
  healthy = rbind(healthy, arrow::read_parquet(filelist[i],
                                               col_select = c("henkilotunnus", "time_to_dg", "tutkimus_lyhenne_yksikko", "tulos", "viitearvot")) %>%
    dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "L-|fP|S-|Pt|Retic|Crea|B-La|P-")) %>%
      mutate(viitearvo_min = gsub("\\/[[:print:]]*", "", gsub("-[[:print:]]*", "", viitearvot)),
             viitearvo_min = ifelse(str_detect(viitearvo_min, "<"), 0, viitearvo_min),
             viitearvo_max = gsub("\\/[[:print:]]*", "", gsub("[[:print:]]{1,3}-", "", viitearvot)),
             viitearvo_max = gsub("<", "", viitearvo_max)))
}

# Load test data
disease = c("de_novo_AML", "MDS", "MF")
results1 = data.frame()

# Loop over diseases
for (i in disease) {
  
  # Load test data
  df = read.csv(paste0("mounts/research/husdatalake/disease/scripts/Preleukemia/oona_new/trajectory_model/XGBoost_COX/youden/", i, "_test_data_with_predictions.csv")) %>%
    dplyr::select(henkilotunnus=patient_id, time_to_dg, disease = disease_status, predicted_label) %>%
    mutate(time_to_dg = ifelse(disease == 1, -time_to_dg, time_to_dg)) %>%
    distinct()
  
  # Real lab values
  ## Patients
  lab = readRDS(paste0(export, "/", gsub("_", " ", i), "_lab_demo.rds")) %>%
    dplyr::select(henkilotunnus, time_to_dg, tutkimus_lyhenne_yksikko, tulos, viitearvot) %>%
    dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "L-|fP|S-|Pt|Retic|Crea|B-La|P-")) %>%
    mutate(viitearvo_min = gsub("\\/[[:print:]]*", "", gsub("-[[:print:]]*", "", viitearvot)),
           viitearvo_min = ifelse(str_detect(viitearvo_min, "<"), 0, viitearvo_min),
           viitearvo_min = ifelse(is.na(viitearvo_min) & tutkimus_lyhenne_yksikko=="B-Erblast (E9/l)", 0, viitearvo_min),
           viitearvo_max = gsub("\\/[[:print:]]*", "", gsub("[[:print:]]{1,3}-", "", viitearvot)),
           viitearvo_max = gsub("<", "", viitearvo_max),
           viitearvo_max = ifelse(is.na(viitearvo_max) & tutkimus_lyhenne_yksikko=="B-Erblast (E9/l)", 0.01, viitearvo_max))
  
  # Join
  df = df %>%
    dplyr::left_join(bind_rows(lab, healthy)) %>%
    dplyr::filter(!is.na(tulos)) %>%
    distinct() %>%
    dplyr::mutate(TP = ifelse(disease == 1 & predicted_label == 1, 1, 0),
                  FN = ifelse(disease == 1 & predicted_label == 0, 1, 0),
                  FP = ifelse(disease == 0 & predicted_label == 1, 1, 0))
  
  # Lab values
  df = df %>%
    mutate(tulos = tulos-0.01)
  df = df %>%
    mutate(out = ifelse(tulos > viitearvo_max | tulos < viitearvo_min, 1, 0))
  
  # Summarise
  tt = df %>%
    group_by(disease, TP, FN, FP) %>%
    summarise(n = sum(out, na.rm = TRUE),
              all=n(),
              prop = round(100*n/all, 2)) %>%
    ungroup() %>%
    mutate(dg = i)
  tt1 = df %>%
    group_by(disease) %>%
    summarise(n = sum(out, na.rm = TRUE),
              all=n(),
              prop = round(100*n/all, 2)) %>%
    ungroup() %>%
    dplyr::filter(!disease == 0) %>%
    mutate(dg = i)
  tt2 = df %>%
    group_by(disease) %>%
    summarise(n = sum(out, na.rm = TRUE),
              all=n(),
              prop = round(100*n/all, 2)) %>%
    ungroup() %>%
    dplyr::filter(disease == 0) %>%
    mutate(dg = i)
  
  # Join
  tt = bind_rows(tt, tt1)
  tt = bind_rows(tt, tt2)
  
  # Collect
  results1 = rbind(results1, tt)
  
}

# Save
write_xlsx(results1, paste0(results, "/lab_value_reference_range.xlsx"))
