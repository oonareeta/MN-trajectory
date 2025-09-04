# Comorbidities of healthy subjects

# Libraries
source("mounts/research/src/Rfunctions/library.R")


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")


# Read data of controls
df_healthy = data.frame()
for (i in 1:100) {
  print(i)
  df_healthy1 = readRDS(paste0(export, "/healthy/lagged_data_", i, "_.rds")) %>%
    dplyr::select(henkilotunnus, naytteenotto_hetki, time_to_dg) %>%
    distinct()
  df_healthy = rbind(df_healthy1, df_healthy)
}
saveRDS(df_healthy, paste0(export, "/healthy_lab_dates.rds"))


# Filter only healthy patients used in modelling
lagged_data1 = readRDS(paste0(export, "/healthy/lagged_data_final.rds")) %>%
  dplyr::select(henkilotunnus, time_to_dg)
lagged_data1$time_to_dg = -lagged_data1$time_to_dg
df_healthy = df_healthy %>%
  dplyr::inner_join(lagged_data1) %>%
  distinct()
rm(lagged_data1)

# Add mortality
asiakas = arrow::read_parquet("mounts/research/husdatalake/data/demographics/cressida_demographics.parquet")
df_healthy = df_healthy %>%
  dplyr::left_join(asiakas %>%
                     dplyr::select(henkilotunnus, kuolinaika_pvm))
tt = df_healthy %>%
  dplyr::mutate(death = ifelse(is.na(kuolinaika_pvm), 0, 1)) %>%
  distinct(henkilotunnus, death)
table(tt$death); table(tt$death) / nrow(tt)

# Never hospitalized
tt = df_healthy %>%
  dplyr::filter(!henkilotunnus %in% inpatient$henkilotunnus) %>%
  distinct(henkilotunnus)
print(paste0("There are ", nrow(tt), " control subjects who have never been hospitalized"))
df_healthy = df_healthy %>%
  dplyr::filter(!henkilotunnus %in% tt$henkilotunnus)


# Collect hospital location data
dir.create(paste0(export, "/healthy/location/"))
listfiles = list.files("mounts/research/husdatalake/data/multilab/main", full.names = TRUE)
for (i in 1:length(listfiles)) {
  print(i)
  tt = arrow::read_parquet(listfiles[i]) %>%
    dplyr::distinct(henkilotunnus, naytteenotto_hetki, pyytava_toimipiste_selite)
  arrow::write_parquet(tt, paste0(export, "/healthy/location/healthy_lab_location_", i, ".parquet"))
}
## Combine
tt = data.frame()
listfiles = list.files("mounts/research/husdatalake/disease/processed_data/Preleukemia/healthy/location", full.names = TRUE)
for (i in 1:length(listfiles)) {
  print(i)
  tt1 = arrow::read_parquet(listfiles[i]) %>%
    dplyr::distinct(henkilotunnus, naytteenotto_hetki, pyytava_toimipiste_selite)
  tt = bind_rows(tt, tt1)
}

# Filter study controls
tt = tt %>%
  dplyr::filter(henkilotunnus %in% unique(lagged_data1$henkilotunnus))

# Save
arrow::write_parquet(tt, paste0(export, "/healthy/healthy_lab_location.parquet"))
unlink(paste0(export, "/healthy/location"))


# Read data
tt = arrow::read_parquet(paste0(export, "/healthy/healthy_lab_location.parquet"))

print(paste0("Number of controls is ", length(unique(tt$henkilotunnus))))

# Annotate hospitalization
tt = tt %>%
  distinct() %>%
  dplyr::mutate(hospitalization = ifelse(str_detect(pyytava_toimipiste_selite, "(?i) OS|OSAS"), 1, 0))
ta = as.data.frame(table(tt$pyytava_toimipiste_selite))
table(tt$hospitalization)
table(tt$hospitalization)[2]/(nrow(tt))

