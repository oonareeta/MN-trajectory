# Comorbidities of healthy subjects

# Libraries
source("mounts/research/src/Rfunctions/library.R")


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Read data of controls
lagged_data1 = fread(paste0(export, "/lagged_data_MF3.csv"), select = c("henkilotunnus", "disease")) %>%
  distinct(henkilotunnus, .keep_all=TRUE) %>%
  dplyr::filter(disease == 0)

# Collect their diagnoses
if (exists("dg")) { rm(dg) }
filenames = list.files("mounts/research/husdatalake/data/diagnoses/diagnoses", full.names = TRUE, recursive = TRUE, all.files = TRUE)
for (i in 1:length(filenames)) {
  print(paste0(i, "/", length(filenames)))
  dg1 <- arrow::read_parquet(filenames[i], col_select = c(henkilotunnus, diagnoosi_koodi, diagnoosi_selite, kontakti_alku_hetki)) %>%
    dplyr::filter(henkilotunnus %in% lagged_data1$henkilotunnus) %>%
    dplyr::arrange(kontakti_alku_hetki) %>%
    group_by(henkilotunnus, diagnoosi_koodi) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    distinct(henkilotunnus, diagnoosi_koodi, diagnoosi_selite, kontakti_alku_hetki)
  if (i == 1) {
    dg = dg1
  } else {
    dg = rbind(dg, dg1)
  }
}

# Save
saveRDS(dg, paste0(export, "/diagnoses_healthy.rds"))
demo = readRDS(paste0(export, "/demo_healthy.rds"))

# Collect unique diagnoses
dgs1 = dg %>%
  group_by(henkilotunnus, diagnoosi_koodi) %>%
  slice(1) %>%
  ungroup()


# Diagnoses in general
tmp = as.data.frame(table(dgs1$diagnoosi_koodi)) %>%
  dplyr::rename(diagnoosi_koodi = Var1) %>%
  dplyr::left_join(dg %>%
                     dplyr::select(diagnoosi_koodi, diagnoosi_selite) %>%
                     group_by(diagnoosi_koodi) %>%
                     dplyr::slice(1) %>%
                     ungroup())

# Number of diagnoses per patients
## Patients with a hematological condition
tt = dgs1 %>%
  dplyr::filter(str_detect(diagnoosi_koodi, "D(5|6|7|8|9)")) %>%
  distinct(henkilotunnus) %>%
  nrow()
print(paste0(length(unique(dgs1$henkilotunnus)) + nrow(pts_wo_dgs) - tt, " patients with a hematological condition"))
tt = dgs1 %>%
  dplyr::filter(str_detect(diagnoosi_koodi, "D(5|6|7|8|9)")) %>%
  distinct(henkilotunnus, .keep_all = TRUE)
tt = as.data.frame(table(tt$diagnoosi_selite))
## Patients wo any diagnoses
pts_wo_dgs = lagged_data1 %>%
  dplyr::filter(!henkilotunnus %in% dgs1$henkilotunnus) %>%
  dplyr::select(-disease) %>%
  mutate(n = 0)
dg_per_pt = dgs1 %>%
  mutate(diagnoosi_koodi = gsub("\\.[[:print:]]*", "", diagnoosi_koodi)) %>%
  dplyr::filter(nchar(diagnoosi_koodi) < 4) %>%
  distinct(henkilotunnus, diagnoosi_koodi) %>%
  dplyr::group_by(henkilotunnus) %>%
  summarise(n = n()) %>%
  ungroup() %>%
  dplyr::bind_rows(pts_wo_dgs)
print(paste0("Diagnose frequencies")); quantile(dg_per_pt$n)
print(paste0("There are ", nrow(pts_wo_dgs), " patients with 0 diagnoses."))


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
  dplyr::filter(str_detect(diagnoosi_selite, "hyytym|Willebrand")) %>%
  distinct(henkilotunnus, .keep_all = TRUE)
unique(coa$diagnoosi_selite)
print(paste0("There are ", nrow(coa), " patients with coagulation disorder"))
neu = dgs1 %>%
  dplyr::filter(str_detect(diagnoosi_selite, "neut")) %>%
  distinct(henkilotunnus, .keep_all = TRUE)
unique(neu$diagnoosi_selite)
print(paste0("There are ", nrow(neu), " patients with neutrophil disorder"))

length(unique(dgs1$henkilotunnus))

aa = as.data.frame(table(dgs1$diagnoosi_selite))

