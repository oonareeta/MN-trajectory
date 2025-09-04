# Calculate missing data

# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

tt1 = data.frame()

# Loop
for (i in c("MF", "MDS", "de novo AML")) {
  
  # Load data
  lagged_data1 = readRDS(paste0(export, "/lagged_data_", i, "2.rds"))
  
  # Remove data 5-10 years before diagnosis
  lagged_data1 <- lagged_data1 %>%
    dplyr::filter(disease == 0 | time_to_dg < 1824)
  
  ### Replace NaN with NA
  lagged_data1 <- lagged_data1 %>%
    dplyr::select(henkilotunnus, time_to_dg, disease, event_1y, sukupuoli_selite, age, everything())
  lagged_data1 <- lagged_data1 %>%
    mutate_all(~ifelse(is.nan(.), NA, .))
  
  # Variables
  vars = names(lagged_data1)[grep(pattern = "tulos_norm", x = names(lagged_data1))]
  vars = vars[grep(pattern = "l_blast|atyp|rdw_sd|ind|p_kol|trfesat|p_fe", x = vars, invert = TRUE)]
  
  # Cound missing rows
  tt = as.data.frame(sapply(vars, function(x) 100-(round(100/nrow(lagged_data1)*nrow(lagged_data1[!is.na(lagged_data1[[x]]),]),1))))
  colnames(tt)[1] = "Missing"
  tt$Disease = i
  
  # Join
  tt1 = rbind(tt1, tt)
}

# Rowname to column
tt1 = tt1 %>%
  rownames_to_column(var = "Labtest") %>%
  mutate(Labtest = gsub("norm[[:digit:]]", "norm", Labtest))

# Dcast
tt2 = reshape2::dcast(tt1, Labtest~Disease, value.var = "Missing")

# Export
writexl::write_xlsx(tt2, paste0(export, "/missing_data.xlsx"))
