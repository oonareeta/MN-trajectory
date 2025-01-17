# Comorbidities of healthy subjects
# Author: Oscar Brück

# Libraries
source("./library.R")


# General parameters
source("./parameters")


############################################################


# Read data of controls
lagged_data1 = readRDS(paste0(export, "/lagged_data_MF3.rds")) %>%
  distinct(henkilotunnus, .keep_all=TRUE) %>%
  dplyr::filter(disease == 0)

# Collect their diagnoses
dgs = data.frame()
for (dg1 in list.files("mounts/research/husdatalake/data/diagnoses/diagnoses/rds", recursive = TRUE, full.names = TRUE)) {
  print(dg1)
  dgs = rbind(dgs, readRDS(dg1) %>%
                dplyr::filter(henkilotunnus %in% lagged_data1$henkilotunnus) %>%
                dplyr::distinct(henkilotunnus, diagnoosi_koodi, diagnoosi_selite, paadiagnoosi, kontakti_alku_hetki))
}

# Save
saveRDS(dgs, paste0(export, "/diagnoses_healthy.rds"))

# Collect unique diagnoses
dgs1 = dgs %>%
  group_by(henkilotunnus, diagnoosi_koodi) %>%
  slice(1) %>%
  ungroup()


# Look for hematological diagnoses
anemia = dgs1 %>%
  dplyr::filter(str_detect(diagnoosi_selite, "(?i)anemia")) %>%
  distinct(henkilotunnus)
print(paste0("There are ", nrow(anemia), " patients with anemia"))
plt = dgs1 %>%
  dplyr::filter(str_detect(diagnoosi_selite, "trombosyto|verihiuta")) %>%
  distinct(henkilotunnus, .keep_all = TRUE)
unique(plt$diagnoosi_selite)
print(paste0("There are ", nrow(plt), " patients with platelet disorders"))
coa = dgs1 %>%
  dplyr::filter(str_detect(diagnoosi_selite, "hyytym")) %>%
  distinct(henkilotunnus, .keep_all = TRUE)
unique(coa$diagnoosi_selite)
print(paste0("There are ", nrow(coa), " patients with coagulation disorder"))
neu = dgs1 %>%
  dplyr::filter(str_detect(diagnoosi_selite, "neut")) %>%
  distinct(henkilotunnus, .keep_all = TRUE)
unique(neu$diagnoosi_selite)
print(paste0("There are ", nrow(neu), " patients with neutrophil disorder"))

tmp = as.data.frame(table(dgs2$diagnoosi_koodi, dgs2$diagnoosi_selite)) %>%
  dplyr::filter(Freq > 0)
length(unique(dgs1$henkilotunnus))

aa = as.data.frame(table(dgs1$diagnoosi_selite))
