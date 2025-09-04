# Measure accuracy of imputation

# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")


# Read data before imputation
i = "MDS"
lagged_data1 = readRDS(paste0(export, "/lagged_data_", i, "2.rds"))
i = "MF"
lagged_data1 = rbind(lagged_data1, readRDS(paste0(export, "/lagged_data_", i, "2.rds")))
i = "de novo AML"
## Join
lagged_data1 = rbind(lagged_data1, readRDS(paste0(export, "/lagged_data_", i, "2.rds")))
lagged_data1$disease = ifelse(lagged_data1$disease %in% c("MDS", "MF", "de novo AML"), 1, 0)
## Save
lagged_data2 = lagged_data1

# Remove redundant variables
lagged_data2 = lagged_data2 %>%
  dplyr::mutate(time_to_dg = -time_to_dg)

# Calculate how many times lab tests have been taken within the last month
lagged_data2 <- lagged_data2 %>%
  arrange(henkilotunnus, desc(time_to_dg)) %>%
  group_by(henkilotunnus) %>%
  mutate(rows_in_last_month = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 31 & time_to_dg >= t)),
         rows_in_last_year = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 366 & time_to_dg >= t))) %>%
  ungroup()

# Remove data 5-10 years before diagnosis
lagged_data2 <- lagged_data2 %>%
  dplyr::filter(disease == 0 | time_to_dg < 1824)

# Remove atyply and i_ind if still there
lagged_data2 = lagged_data2 %>%
  dplyr::select(-all_of(names(lagged_data2)[grep(pattern = "i_ind", ignore.case = TRUE, names(lagged_data2))])) %>%
  dplyr::select(-all_of(names(lagged_data2)[grep(pattern = "l_atyp_ly_percent", ignore.case = TRUE, names(lagged_data2))]))

# Manually impute NA to 0.01 for erblast as these are by default 0 if not reported
lagged_data2$b_erblast_e9_l_tulos_norm = ifelse(is.na(lagged_data2$b_erblast_e9_l_tulos_norm), 0.01, lagged_data2$b_erblast_e9_l_tulos_norm)
lagged_data2$b_erblast_e9_l_mean_last_values_365_days = ifelse(is.na(lagged_data2$b_erblast_e9_l_mean_last_values_365_days), 0.01, lagged_data2$b_erblast_e9_l_mean_last_values_365_days)
lagged_data2$b_erblast_e9_l_mean_last_values_1825_days = ifelse(is.na(lagged_data2$b_erblast_e9_l_mean_last_values_1825_days), 0.01, lagged_data2$b_erblast_e9_l_mean_last_values_1825_days)
lagged_data2$b_erblast_e9_l_mean_last_values_1095_days = ifelse(is.na(lagged_data2$b_erblast_e9_l_mean_last_values_1095_days), 0.01, lagged_data2$b_erblast_e9_l_mean_last_values_1095_days)
lagged_data2$b_erblast_e9_l_trend_from_365_d = ifelse(is.na(lagged_data2$b_erblast_e9_l_trend_from_365_d), 0, lagged_data2$b_erblast_e9_l_trend_from_365_d)
lagged_data2$b_erblast_e9_l_trend_from_1095_d = ifelse(is.na(lagged_data2$b_erblast_e9_l_trend_from_1095_d), 0, lagged_data2$b_erblast_e9_l_trend_from_1095_d)
lagged_data2$b_erblast_e9_l_trend_from_1825_d = ifelse(is.na(lagged_data2$b_erblast_e9_l_trend_from_1825_d), 0, lagged_data2$b_erblast_e9_l_trend_from_1825_d)

## Remove if fewer than 3 rows per patient
lagged_data1_sr = lagged_data2 %>%
  group_by(henkilotunnus) %>%
  summarise(n = n()) %>%
  dplyr::filter(n >= 3)
lagged_data3 = lagged_data2 %>%
  dplyr::filter(henkilotunnus %in% lagged_data1_sr$henkilotunnus)




# Repeat for healthy
lagged_data1= data.frame()
for (i in 1:10) {
  lagged_data1 = bind_rows(lagged_data1, readRDS(paste0(export, "/healthy/lagged_data_", i, "_2.rds")))
}
lagged_data1$disease = 0

# Remove redundant variables
lagged_data1 = lagged_data1 %>%
  dplyr::mutate(time_to_dg = -time_to_dg)

# Calculate how many times lab tests have been taken within the last month
lagged_data1 <- lagged_data1 %>%
  # slice(1:20000) %>%
  arrange(henkilotunnus, desc(time_to_dg)) %>%
  group_by(henkilotunnus) %>%
  mutate(rows_in_last_month = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 31 & time_to_dg >= t)),
         rows_in_last_year = sapply(time_to_dg, function(t) sum(time_to_dg <= t + 366 & time_to_dg >= t))) %>%
  ungroup()

# Remove data 5-10 years before diagnosis
lagged_data1 <- lagged_data1 %>%
  dplyr::filter(disease == 0 | time_to_dg < 1824)

# Remove atyply and i_ind if still there
lagged_data1 = lagged_data1 %>%
  dplyr::select(-all_of(names(lagged_data1)[grep(pattern = "i_ind", ignore.case = TRUE, names(lagged_data1))])) %>%
  dplyr::select(-all_of(names(lagged_data1)[grep(pattern = "l_atyp_ly_percent", ignore.case = TRUE, names(lagged_data1))]))

# Manually impute NA to 0.01 for erblast as these are by default 0 if not reported
lagged_data1$b_erblast_e9_l_tulos_norm = ifelse(is.na(lagged_data1$b_erblast_e9_l_tulos_norm), 0.01, lagged_data1$b_erblast_e9_l_tulos_norm)
lagged_data1$b_erblast_e9_l_mean_last_values_365_days = ifelse(is.na(lagged_data1$b_erblast_e9_l_mean_last_values_365_days), 0.01, lagged_data1$b_erblast_e9_l_mean_last_values_365_days)
lagged_data1$b_erblast_e9_l_mean_last_values_1825_days = ifelse(is.na(lagged_data1$b_erblast_e9_l_mean_last_values_1825_days), 0.01, lagged_data1$b_erblast_e9_l_mean_last_values_1825_days)
lagged_data1$b_erblast_e9_l_mean_last_values_1095_days = ifelse(is.na(lagged_data1$b_erblast_e9_l_mean_last_values_1095_days), 0.01, lagged_data1$b_erblast_e9_l_mean_last_values_1095_days)
lagged_data1$b_erblast_e9_l_trend_from_365_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_365_d), 0, lagged_data1$b_erblast_e9_l_trend_from_365_d)
lagged_data1$b_erblast_e9_l_trend_from_1095_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_1095_d), 0, lagged_data1$b_erblast_e9_l_trend_from_1095_d)
lagged_data1$b_erblast_e9_l_trend_from_1825_d = ifelse(is.na(lagged_data1$b_erblast_e9_l_trend_from_1825_d), 0, lagged_data1$b_erblast_e9_l_trend_from_1825_d)



# Combine healthy and disease
lagged_data4 = bind_rows(lagged_data1, lagged_data3)

# Select henkilotunnus and lab values
lagged_data4 = lagged_data4 %>%
  dplyr::select(1, 2, 9:ncol(.)) %>%
  dplyr::filter(time_to_dg >= 90)
# Melt
lagged_data4 = lagged_data4 %>%
  reshape2::melt(id.vars = c("henkilotunnus", "time_to_dg"))
# Remove NA
lagged_data4 = lagged_data4 %>%
  dplyr::filter(!is.na(value))
lagged_data4 = lagged_data4 %>%
  dplyr::filter(!value=="NaN")
# Remove val variables
lagged_data4 = lagged_data4 %>%
  dplyr::filter(!str_detect(variable, "_val_"))
# Export
fwrite(lagged_data4, paste0(export, "/impute_before.csv"))
lagged_data4 = fread(paste0(export, "/impute_before.csv"))


# Replace 10% of 'value' with NA using dplyr
lagged_data4 <- lagged_data4 %>%
  distinct() %>%
  mutate(value_impute = if_else(row_number() %in% sample(n(), size = floor(0.2 * n())), NA_real_, value))

# Impute
lagged_data4 <- lagged_data4 %>%
  dplyr::arrange(desc(time_to_dg)) %>%
  group_by(henkilotunnus, variable) %>%
  
  # Imputation using interpolation
  mutate(value_impute1 = zoo::na.approx(value_impute, na.rm = FALSE))

# Statistics
cor.test(lagged_data4[is.na(lagged_data4$value_impute),]$value_impute1, lagged_data4[is.na(lagged_data4$value_impute),]$value, method = "spearman")

# Export
fwrite(lagged_data4, paste0(export, "/impute_after.csv"))
lagged_data4 = fread(paste0(export, "/impute_after.csv"))


# Plot
tt = lagged_data4 %>%
  dplyr::filter(is.na(value_impute))
cor.test(tt$value_impute1, tt$value, method = "spearman")

g = ggplot(tt, aes(x = value_impute1, y = value)) +
  geom_point(size = 0.2, color = "black") +
  geom_smooth(method = "lm") +
  labs(y="Imputed normalized laboratory value",
       x="Original normalized laboratory value") +
  stat_cor(method = "spearman",
           label.x = -1000,
           label.y = 2000) +
  xlim(-1000, 2000) +
  ylim(-1000, 2000) +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_text(size=14, colour = "black"),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=14, face="bold", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"); g

# Export
ggsave(plot = g,
       filename = paste0(results, "/Imputation_quality.png"),
       height = 4, width = 4, units = "in", dpi = 300)
