# Heatmap

# Libraries
source("mounts/research/src/Rfunctions/library.R")


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")


##### Labs for MN patients with sets #####


diseases = c("de novo AML", "MDS", "MF") 
tt1_sr = data.frame()

# Collect all information at once
for (i in diseases) {
  
  print(i)
  
  
  # Read data
  x_train = fread(paste0(export, "/lab_data_for_modelling_", i, ".csv")) %>%
    dplyr::filter(time_to_dg < -90)
  sets = fread(paste0(export, "/cohort_ids/", gsub(" ", "_", i), "_cohort_ids.csv"))
  
  
  # Select patients included in the model training or evaluation
  x_train = x_train %>%
    dplyr::inner_join(sets %>%
                        dplyr::rename(henkilotunnus = pt_id))
  
  # Remove unnecessary variables
  x_train = x_train %>%
    dplyr::select(henkilotunnus, time_to_dg, tutkimus_lyhenne_yksikko, tulos, sukupuoli_selite, age, set)
  
  
  # Summarise
  tt1 = x_train %>%
    mutate(tulos = tulos-0.01) %>%
    group_by(tutkimus_lyhenne_yksikko, set) %>%
    summarise(tulos_m = median(tulos, na.rm=TRUE),
              tulos_q1 = quantile(tulos, na.rm=TRUE, probs=0.25),
              tulos_q3 = quantile(tulos, na.rm=TRUE, probs=0.75)) %>%
    ungroup() %>%
    mutate(tulos = paste0(tulos_m, " (", tulos_q1, "-", tulos_q3, ")"),
           disease = i) %>%
    dplyr::select(tutkimus_lyhenne_yksikko, disease, set, tulos)
  
  # Join
  tt1_sr = rbind(tt1_sr, tt1) %>%
    dplyr::arrange(set, disease)
  
}

# Dcast
tt1_sr1 = tt1_sr %>%
  distinct() %>%
  reshape2::dcast(tutkimus_lyhenne_yksikko + set ~ disease, value.var = "tulos") %>%
  dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "^S-|^Pt-|^P-|Segm|Myel|Metam|L-Eryt|L-Blas|Band|^fP|Retic|^B-La"))

# Show results
view(tt1_sr1)

# Export
fwrite(tt1_sr1, paste0(export, "/Patient_lab_table_MN.csv"))


##### Labs for MN patients without sets #####


diseases = c("de novo AML", "MDS", "MF") 
tt1_sr = data.frame()

# Collect all information at once
for (i in diseases) {
  
  print(i)
  
  
  # Read data
  x_train = fread(paste0(export, "/lab_data_for_modelling_", i, ".csv")) %>%
    dplyr::filter(time_to_dg < -90)
  sets = fread(paste0(export, "/cohort_ids/", gsub(" ", "_", i), "_cohort_ids.csv"))
  
  
  # Select patients included in the model training or evaluation
  x_train = x_train %>%
    dplyr::filter(henkilotunnus %in% sets$pt_id)
  
  # Remove unnecessary variables
  x_train = x_train %>%
    dplyr::select(henkilotunnus, time_to_dg, tutkimus_lyhenne_yksikko, tulos, sukupuoli_selite, age)
  
  
  # Summarise
  tt1 = x_train %>%
    mutate(tulos = tulos-0.01) %>%
    group_by(tutkimus_lyhenne_yksikko) %>%
    summarise(tulos_m = median(tulos, na.rm=TRUE),
              tulos_q1 = quantile(tulos, na.rm=TRUE, probs=0.25),
              tulos_q3 = quantile(tulos, na.rm=TRUE, probs=0.75)) %>%
    ungroup() %>%
    mutate(tulos = paste0(tulos_m, " (", tulos_q1, "-", tulos_q3, ")"),
           disease = i) %>%
    dplyr::select(tutkimus_lyhenne_yksikko, disease, tulos)
  
  # Join
  tt1_sr = rbind(tt1_sr, tt1) %>%
    dplyr::arrange(disease)
  
}

# Dcast
tt1_sr1 = tt1_sr %>%
  distinct() %>%
  reshape2::dcast(tutkimus_lyhenne_yksikko ~ disease, value.var = "tulos") %>%
  dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "^S-|^Pt-|^P-|Segm|Myel|Metam|L-Eryt|L-Blas|Band|^fP|Retic|^B-La"))

# Show results
view(tt1_sr1)

# Export
fwrite(tt1_sr1, paste0(export, "/Patient_lab_table_MN_total.csv"))


##### Labs for controls with sets #####


x_train = data.frame()
filelist1 = list.files(paste0(export, "/healthy/"), pattern = "lab_data_for_modelling_")
# Collect all information at once
for (i in 1:length(filelist1)) {
  
  print(i)
  
  # Read data
  x_train1 = fread(paste0(export, "/healthy/lab_data_for_modelling_", i, ".csv"),
                   select = c("henkilotunnus", "time_to_dg", "tutkimus_lyhenne_yksikko", "tulos", "sukupuoli_selite", "age")) %>%
    dplyr::filter(time_to_dg <= -365)
  
  # Join
  x_train = rbind(x_train, x_train1)
  
}

# Read sets
sets = fread(paste0(export, "/cohort_ids/MDS_cohort_ids.csv"))


# Select patients included in the model training or evaluation
x_train = x_train %>%
  dplyr::inner_join(sets %>%
                      dplyr::rename(henkilotunnus = pt_id))


# Summarise
tt1 = x_train %>%
  mutate(tulos = tulos-0.01,
         tutkimus_lyhenne_yksikko = gsub("\\/solu", "", 
                                         gsub("PG", "pg", 
                                              gsub("FL", "fl", 
                                                   gsub("Mg", "mg", 
                                                        gsub("09\\/", "9/", 
                                                             gsub("MM\\/", "mm/", 
                                                                  gsub("\\/H", "/h", 
                                                                       gsub("\\/l", "/L", 
                                                                            gsub("G\\/", "g/", tutkimus_lyhenne_yksikko)))))))))) %>%
  group_by(tutkimus_lyhenne_yksikko, set) %>%
  summarise(tulos_m = median(tulos, na.rm=TRUE),
            tulos_q1 = quantile(tulos, na.rm=TRUE, probs=0.25),
            tulos_q3 = quantile(tulos, na.rm=TRUE, probs=0.75)) %>%
  ungroup() %>%
  mutate(tulos = paste0(tulos_m, " (", tulos_q1, "-", tulos_q3, ")"),
         disease = "Control") %>%
  dplyr::select(tutkimus_lyhenne_yksikko, disease, set, tulos)
# Join
tt1_sr = tt1 %>%
  dplyr::arrange(set)

# Dcast
tt1_sr1 = tt1_sr %>%
  distinct() %>%
  reshape2::dcast(tutkimus_lyhenne_yksikko + set ~ disease, value.var = "tulos") %>%
  dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "^S-|^Pt-|^P-|Segm|Myel|Metam|L-Eryt|L-Blas|Band|^fP|Retic|^B-La"))

# Show results
view(tt1_sr1)


# Export
fwrite(tt1_sr1, paste0(export, "/Patient_lab_table_Controls.csv"))


##### Labs for controls without sets #####


x_train = data.frame()
filelist1 = list.files(paste0(export, "/healthy/"), pattern = "lab_data_for_modelling_")
# Collect all information at once
for (i in 1:length(filelist1)) {
  
  print(i)
  
  # Read data
  x_train1 = fread(paste0(export, "/healthy/lab_data_for_modelling_", i, ".csv"),
                   select = c("henkilotunnus", "time_to_dg", "tutkimus_lyhenne_yksikko", "tulos", "sukupuoli_selite", "age")) %>%
    dplyr::filter(time_to_dg <= -365)
  
  # Join
  x_train = rbind(x_train, x_train1)
  
}

# Read sets
sets = fread(paste0(export, "/cohort_ids/MDS_cohort_ids.csv"))


# Select patients included in the model training or evaluation
x_train = x_train %>%
  dplyr::filter(henkilotunnus %in% sets$pt_id)


# Summarise
tt1 = x_train %>%
  mutate(tulos = tulos-0.01,
         tutkimus_lyhenne_yksikko = gsub("\\/solu", "", 
                                         gsub("PG", "pg", 
                                              gsub("FL", "fl", 
                                                   gsub("Mg", "mg", 
                                                        gsub("09\\/", "9/", 
                                                             gsub("MM\\/", "mm/", 
                                                                  gsub("\\/H", "/h", 
                                                                       gsub("\\/l", "/L", 
                                                                            gsub("G\\/", "g/", tutkimus_lyhenne_yksikko)))))))))) %>%
  group_by(tutkimus_lyhenne_yksikko) %>%
  summarise(tulos_m = median(tulos, na.rm=TRUE),
            tulos_q1 = quantile(tulos, na.rm=TRUE, probs=0.25),
            tulos_q3 = quantile(tulos, na.rm=TRUE, probs=0.75)) %>%
  ungroup() %>%
  mutate(tulos = paste0(tulos_m, " (", tulos_q1, "-", tulos_q3, ")"),
         disease = "Control") %>%
  dplyr::select(tutkimus_lyhenne_yksikko, disease, tulos)

# Join
tt1_sr1 = tt1

# Dcast
tt1_sr1 = tt1_sr1 %>%
  distinct() %>%
  reshape2::dcast(tutkimus_lyhenne_yksikko ~ disease, value.var = "tulos") %>%
  dplyr::filter(!str_detect(tutkimus_lyhenne_yksikko, "^S-|^Pt-|^P-|Segm|Myel|Metam|L-Eryt|L-Blas|Band|^fP|Retic|^B-La"))

# Show results
view(tt1_sr1)


# Export
fwrite(tt1_sr1, paste0(export, "/Patient_lab_table_Controls_Total.csv"))


############### General information by sets ##################################


# Disease
tt2_sr = data.frame()
tt3_sr = data.frame()
for (i in diseases) {
  
  print(i)
  
  rm(x1)
  
  # Time to event
  x1 = fread(paste0("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_", i, "3.csv")) %>%
    dplyr::filter(disease == 1)
  # Select patients included in the model training or evaluation
  sets = fread(paste0(export, "/cohort_ids/", gsub(" ", "_", i), "_cohort_ids.csv"))
  x1 = x1 %>%
    dplyr::filter(time_to_dg > 90) %>%
    dplyr::inner_join(sets %>%
                        dplyr::rename(henkilotunnus = pt_id)) %>%
    mutate(time_to_dg = -time_to_dg)
  
  
  # Summarise
  tt2 = x1 %>%
    arrange(time_to_dg) %>%
    group_by(henkilotunnus) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    group_by(set) %>%
    summarise(time_to_dg_m = round(median(time_to_dg)/365.24, 2),
              time_to_dg_q1 = round(quantile(time_to_dg, probs=0.25)/365.24, 2),
              time_to_dg_q3 = round(quantile(time_to_dg, probs=0.75)/365.24, 2),
              male = sum(sukupuoli_selite=="Mies"),
              female = sum(sukupuoli_selite=="Nainen"),
              n = n(),
              male_p = round(100*male/n(), 3),
              female_p = round(100*female/n(), 3),
              age_m = median(age, na.rm=TRUE),
              age_q1 = quantile(age, na.rm=TRUE, probs=0.25),
              age_q3 = quantile(age, na.rm=TRUE, probs=0.75)) %>%
    ungroup() %>%
    mutate(male = paste0(male, "/", n, " (", male_p, ")"),
           female = paste0(female, "/", n, " (", female_p, ")"),
           age = paste0(age_m, " (", age_q1, "-", age_q3, ")"),
           time_to_dg = paste0(time_to_dg_m, " (", time_to_dg_q1, "-", time_to_dg_q3, ")"),
           disease = i) %>%
    dplyr::select(disease, set, age, male, female, n, time_to_dg)
  
  # Summarise
  tt3 = x1 %>%
    distinct(set, henkilotunnus, time_to_dg) %>%
    group_by(set, henkilotunnus) %>%
    summarise(n1 = n()) %>%
    ungroup() %>%
    group_by(set) %>%
    dplyr::summarise(n = median(n1, na.rm=TRUE),
                     n_q1 = quantile(n1, na.rm=TRUE, probs=0.25),
                     n_q3 = quantile(n1, na.rm=TRUE, probs=0.75)) %>%
    ungroup() %>%
    mutate(n = paste0(n, " (", n_q1, "-", n_q3, ")"),
           disease = i) %>%
    dplyr::select(disease, set, n)
  
  # Join
  tt2_sr = rbind(tt2_sr, tt2)
  tt3_sr = rbind(tt3_sr, tt3)
  
}

view(tt2_sr)
view(tt3_sr)


# Controls

# Time to event
x1 = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_de novo AML3.csv") %>%
  dplyr::filter(disease == 0)
# Select patients included in the model training or evaluation
sets = fread(paste0(export, "/cohort_ids/de_novo_AML_cohort_ids.csv"))
x1 = x1 %>%
  dplyr::inner_join(sets %>%
                      dplyr::rename(henkilotunnus = pt_id))


# Summarise
tt2 = x1 %>%
  arrange(time_to_dg) %>%
  group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  group_by(set) %>%
  summarise(time_to_dg_m = round(median(time_to_dg)/365.24, 2),
            time_to_dg_q1 = round(quantile(time_to_dg, probs=0.25)/365.24, 2),
            time_to_dg_q3 = round(quantile(time_to_dg, probs=0.75)/365.24, 2),
            male = sum(sukupuoli_selite=="Mies"),
            female = sum(sukupuoli_selite=="Nainen"),
            n = n(),
            male_p = round(100*male/n(), 3),
            female_p = round(100*female/n(), 3),
            age_m = median(age, na.rm=TRUE),
            age_q1 = quantile(age, na.rm=TRUE, probs=0.25),
            age_q3 = quantile(age, na.rm=TRUE, probs=0.75)) %>%
  ungroup() %>%
  mutate(male = paste0(male, "/", n, " (", male_p, ")"),
         female = paste0(female, "/", n, " (", female_p, ")"),
         age = paste0(age_m, " (", age_q1, "-", age_q3, ")"),
         time_to_dg = paste0(time_to_dg_m, " (", time_to_dg_q1, "-", time_to_dg_q3, ")"),
         disease = "Control") %>%
  dplyr::select(disease, set, age, male, female, n, time_to_dg)

# Summarise
tt3 = x1 %>%
  distinct(set, henkilotunnus, time_to_dg) %>%
  group_by(set, henkilotunnus) %>%
  summarise(n1 = n()) %>%
  ungroup() %>%
  group_by(set) %>%
  dplyr::summarise(n = median(n1, na.rm=TRUE),
                   n_q1 = quantile(n1, na.rm=TRUE, probs=0.25),
                   n_q3 = quantile(n1, na.rm=TRUE, probs=0.75)) %>%
  ungroup() %>%
  mutate(n = paste0(n, " (", n_q1, "-", n_q3, ")"),
         disease = "Control") %>%
  dplyr::select(disease, set, n)

# Save
tt2_controls = tt2
tt3_controls = tt3

view(tt2_controls)
view(tt3_controls)




############### General information without sets ##################################


# Disease
tt2_sr = data.frame()
tt3_sr = data.frame()
for (i in diseases) {
  
  print(i)
  
  rm(x1)
  
  # Time to event
  x1 = fread(paste0("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_", i, "3.csv")) %>%
    dplyr::filter(disease == 1)
  # Select patients included in the model training or evaluation
  x1 = x1 %>%
    dplyr::filter(time_to_dg > 90) %>%
    mutate(time_to_dg = -time_to_dg)
  
  
  # Select patients included in the model training or evaluation
  sets = fread(paste0(export, "/cohort_ids/", gsub(" ", "_", i), "_cohort_ids.csv"))
  x1 = x1 %>%
    dplyr::filter(henkilotunnus %in% sets$pt_id)
  
  
  # Summarise
  tt2 = x1 %>%
    arrange(time_to_dg) %>%
    group_by(henkilotunnus) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    summarise(time_to_dg_m = round(median(time_to_dg)/365.24, 2),
              time_to_dg_q1 = round(quantile(time_to_dg, probs=0.25)/365.24, 2),
              time_to_dg_q3 = round(quantile(time_to_dg, probs=0.75)/365.24, 2),
              male = sum(sukupuoli_selite=="Mies"),
              female = sum(sukupuoli_selite=="Nainen"),
              n = n(),
              male_p = round(100*male/n(), 3),
              female_p = round(100*female/n(), 3),
              age_m = median(age, na.rm=TRUE),
              age_q1 = quantile(age, na.rm=TRUE, probs=0.25),
              age_q3 = quantile(age, na.rm=TRUE, probs=0.75)) %>%
    ungroup() %>%
    mutate(male = paste0(male, "/", n, " (", male_p, ")"),
           female = paste0(female, "/", n, " (", female_p, ")"),
           age = paste0(age_m, " (", age_q1, "-", age_q3, ")"),
           time_to_dg = paste0(time_to_dg_m, " (", time_to_dg_q1, "-", time_to_dg_q3, ")"),
           disease = i) %>%
    dplyr::select(disease, age, male, female, n, time_to_dg)
  
  # Summarise
  tt3 = x1 %>%
    distinct(henkilotunnus, time_to_dg) %>%
    group_by(henkilotunnus) %>%
    summarise(n1 = n()) %>%
    ungroup() %>%
    dplyr::summarise(n = median(n1, na.rm=TRUE),
                     n_q1 = quantile(n1, na.rm=TRUE, probs=0.25),
                     n_q3 = quantile(n1, na.rm=TRUE, probs=0.75)) %>%
    mutate(n = paste0(n, " (", n_q1, "-", n_q3, ")"),
           disease = i) %>%
    dplyr::select(disease, n)
  
  # Join
  tt2_sr = rbind(tt2_sr, tt2)
  tt3_sr = rbind(tt3_sr, tt3)
  
}

view(tt2_sr)
view(tt3_sr)


# Controls

# Time to event
x1 = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_de novo AML3.csv") %>%
  dplyr::filter(disease == 0)


# Select patients included in the model training or evaluation
sets = fread(paste0(export, "/cohort_ids/de_novo_AML_cohort_ids.csv"))
x1 = x1 %>%
  dplyr::filter(henkilotunnus %in% sets$pt_id)


# Summarise
tt2 = x1 %>%
  arrange(time_to_dg) %>%
  group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup() %>%
  summarise(time_to_dg_m = round(median(time_to_dg)/365.24, 2),
            time_to_dg_q1 = round(quantile(time_to_dg, probs=0.25)/365.24, 2),
            time_to_dg_q3 = round(quantile(time_to_dg, probs=0.75)/365.24, 2),
            male = sum(sukupuoli_selite=="Mies"),
            female = sum(sukupuoli_selite=="Nainen"),
            n = n(),
            male_p = round(100*male/n(), 3),
            female_p = round(100*female/n(), 3),
            age_m = median(age, na.rm=TRUE),
            age_q1 = quantile(age, na.rm=TRUE, probs=0.25),
            age_q3 = quantile(age, na.rm=TRUE, probs=0.75)) %>%
  ungroup() %>%
  mutate(male = paste0(male, "/", n, " (", male_p, ")"),
         female = paste0(female, "/", n, " (", female_p, ")"),
         age = paste0(age_m, " (", age_q1, "-", age_q3, ")"),
         time_to_dg = paste0(time_to_dg_m, " (", time_to_dg_q1, "-", time_to_dg_q3, ")"),
         disease = "Control") %>%
  dplyr::select(disease, age, male, female, n, time_to_dg)

# Summarise
tt3 = x1 %>%
  distinct(henkilotunnus, time_to_dg) %>%
  group_by(henkilotunnus) %>%
  summarise(n1 = n()) %>%
  ungroup() %>%
  dplyr::summarise(n = median(n1, na.rm=TRUE),
                   n_q1 = quantile(n1, na.rm=TRUE, probs=0.25),
                   n_q3 = quantile(n1, na.rm=TRUE, probs=0.75)) %>%
  mutate(n = paste0(n, " (", n_q1, "-", n_q3, ")"),
         disease = "Control") %>%
  dplyr::select(disease, n)

# Save
tt2_controls = tt2
tt3_controls = tt3

view(tt2_controls)
view(tt3_controls)

