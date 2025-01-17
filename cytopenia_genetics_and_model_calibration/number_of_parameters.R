# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# de novo AML
x_train = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_de novo AML3.csv") %>%
  dplyr::filter(disease == 0 | (disease==1 & time_to_dg > 30)) %>%
  dplyr::select(-starts_with(c("f_p", "e_retic", "l_blast", "e_rdw_sd")))
# MDS
x_train = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_MDS3.csv") %>%
  dplyr::filter(disease == 0 | (disease==1 & time_to_dg > 30)) %>%
  dplyr::select(-starts_with(c("f_p", "e_retic", "l_blast", "e_rdw_sd")))
# MF
x_train = fread("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_MF3.csv") %>%
  dplyr::filter(disease == 0 | (disease==1 & time_to_dg > 30)) %>%
  dplyr::select(-starts_with(c("f_p", "e_retic", "l_blast", "e_rdw_sd")))

# Remove PV and ET patients identified based on JAK2 mutations and diagnosis table
x_train = x_train %>%
  dplyr::filter(!henkilotunnus %in% c('02139_43396117','02139_84320876','02139_05958098','02139_09649728','02139_43105799','02139_86518942','02139_58286011',
                                      '02139_74536441','02139_04048905','02139_66454572','02139_84083902','02139_08144748','02139_93742663','02139_09325470',
                                      '02139_98606703','02139_32936070','02139_74525166','02139_01261095','02139_72491838','02139_16511870','02139_04981848',
                                      '02139_72056477','02139_47431683','02139_91416073'))

# Unique patients
length(unique(x_train[x_train$disease==0,]$henkilotunnus))
length(unique(x_train[x_train$disease==1,]$henkilotunnus))

# Variable per lab test
names(x_train)[grep(pattern = "monos", x = names(x_train), ignore.case = TRUE)]

# Lab test
names(x_train)[grep(pattern = "mean_last_values_365_days", x = names(x_train), ignore.case = TRUE)]

# Healthy patients
nrow(x_train[x_train$disease==0,])
nrow(x_train[x_train$disease==1,])

# 3 demographic = age, gender, number of labs
# 21 lab test
# 7 variable per lab test
# MF 4506 timepoints
# MDS 12549 timepoints
# de novo AML 5939 timepoints
# healthy 424478 timepoints
(3+21*7)*(5939 + 12549 + 4506 + 424478)

