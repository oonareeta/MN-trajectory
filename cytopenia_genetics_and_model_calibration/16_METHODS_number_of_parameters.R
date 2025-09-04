# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")


# de novo AML
x_train = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_de novo AML3.csv") %>%
  dplyr::filter(disease == 0 | (disease==1 & time_to_dg > 90))
sets = fread(paste0(export, "/cohort_ids/de_novo_AML_cohort_ids.csv"))
# MDS
x_train = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_MDS3.csv") %>%
  dplyr::filter(disease == 0 | (disease==1 & time_to_dg > 90))
sets = fread(paste0(export, "/cohort_ids/MDS_cohort_ids.csv"))
# MF
x_train = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_MF3.csv") %>%
  dplyr::filter(disease == 0 | (disease==1 & time_to_dg > 90))
sets = fread(paste0(export, "/cohort_ids/MF_cohort_ids.csv"))


# Select patients included in the model training or evaluation
x_train = x_train %>%
  dplyr::filter(henkilotunnus %in% sets$pt_id)


# Derivation and test sets
x_train = x_train %>%
  # dplyr::filter(henkilotunnus %in% sets[sets$set=="deriv",]$pt_id) %>%
  # dplyr::filter(henkilotunnus %in% sets[sets$set=="internal test",]$pt_id) %>%
  dplyr::filter(henkilotunnus %in% sets$pt_id) %>%
  filter(disease == 0) %>%
  distinct(henkilotunnus) %>%
  nrow()

# Unique patients
length(unique(x_train[x_train$disease==0,]$henkilotunnus))
length(unique(x_train[x_train$disease==1,]$henkilotunnus))

# Variable per lab test
names(x_train)[grep(pattern = "monos", x = names(x_train), ignore.case = TRUE)]
names(x_train)[grep(pattern = "ferrit", x = names(x_train), ignore.case = TRUE)]

# Lab test
names(x_train)[grep(pattern = "trend_from_365_d", x = names(x_train), ignore.case = TRUE)]
names(x_train)[grep(pattern = "val_365_n", x = names(x_train), ignore.case = TRUE)]

# Healthy patients
nrow(x_train[x_train$disease==0,])
nrow(x_train[x_train$disease==1,])


# 4 demographic = age, gender, number of labs
# 21 lab test with trends x7
# 10 lab test with val_n x3
# 7 variable per lab test
# MF 818 timepoints (n=69)
# MDS 9539 timepoints (n=629)
# de novo AML 4814 timepoints (n=365)
# healthy 4483684 timepoints (n=121739)
(4 + 21*7 + 10*3)*(818 + 9539 + 4814 + 4483675)
31*(818 + 9539 + 4814 + 4483675)

