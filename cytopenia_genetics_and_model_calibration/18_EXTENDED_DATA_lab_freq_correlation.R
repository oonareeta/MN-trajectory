# Calculate laboratory test frequencies

# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Read data
lagged_data1 = fread(paste0(export, "/lagged_data_MDS3.csv"))

# Sample 100 000 rows
lagged_data2 = lagged_data1 %>%
  dplyr::slice_sample(n=100000)

# Plot
g = ggplot(lagged_data2, aes(x = e_mcv_fl_tulos_norm_val_1825_n, y = b_hb_g_l_tulos_norm_val_1825_n)) +
  geom_point(size = 0.2, color = "black") +
  labs(
    y="Freq E-MCV (fl) -5y",
    x="Freq B-Hb (g/L) -5y") +
  stat_cor(method = "spearman") +
  theme_bw() +
  theme(axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title = element_text(size=12, colour = "black"),
        axis.line = element_line(colour = "black"),
        plot.title = element_text(size=12, face="bold", colour = "black"),
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank(),
        legend.position = "none"); g

# Export
ggsave(plot = g,
       filename = paste0(results, "/Lab_freq_correlation.png"),
       height = 4, width = 4, units = "in", dpi = 300)
