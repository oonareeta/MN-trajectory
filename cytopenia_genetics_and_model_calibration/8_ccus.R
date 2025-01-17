# CCUS

# Author: Oscar Brück

# Libraries
source("./library.R")
library(MatchIt)
library(survRM2)

# General parameters
source("./parameters")


################ Load data #################################################################


# Labtests
i = "MDS"
labtests = readxl::read_xlsx(paste0(results, "/", i, "/linear_regression_results_5y.xlsx")) %>%
  distinct(Variable, .keep_all = TRUE) %>%
  dplyr::filter(!str_detect(Variable, "(E-RDW-SD|L-Myelos|L-Band|L-RBC|L-Metam|L-Segment|L-Erytrobl|AtypLy|L-Blast)"))
labtests = c("B-Hb (g/l)", "B-Neut (E9/l)", "B-PLT (E9/l)")#, "L-Neut (%)") #"B-Neut (E9/l)")

# Cytopenic patients
dgs1 = readRDS(paste0(export, "/diagnoses_healthy.rds")) %>%
  group_by(henkilotunnus, diagnoosi_koodi) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::filter(str_detect(diagnoosi_selite, "(?i)(anemia|trombosyto|verihiuta|neut)")) %>%
  distinct(henkilotunnus)

# Disease patients
mf = readRDS(paste0(export, "/lagged_data_MF3.rds"))
mf1 = mf %>%
  dplyr::filter(disease == 1) %>%
  distinct(henkilotunnus)
mf = mf %>%
  dplyr::filter(!henkilotunnus %in% mf1$henkilotunnus)


# Load disease data
df_healthy = fread(paste0(export, "/lab_data_for_modelling_MDS.csv")) %>%
  # Keep only the 11680 healthy controls
  dplyr::filter(henkilotunnus %in% unique(mf$henkilotunnus)) %>%
  # Remove cytopenic controls
  dplyr::filter(tutkimus_lyhenne_yksikko %in% labtests) %>%
  dplyr::mutate(min1 = ifelse((tutkimus_lyhenne_yksikko == "B-Hb (g/l)" & sukupuoli_selite == "Nainen" & tulos < 120) |
                                (tutkimus_lyhenne_yksikko == "B-Hb (g/l)" & sukupuoli_selite == "Mies" & tulos < 130) |
                                (tutkimus_lyhenne_yksikko == "B-PLT (E9/l)" & tulos < 150) |
                                (tutkimus_lyhenne_yksikko == "B-Neut (E9/l)" & tulos < 1.8), 1, 0),
                time_to_dg_mo = ceiling(time_to_dg/365.24)) %>%
  dplyr::arrange(desc(min1)) %>%
  dplyr::mutate(time_to_dg_mo1 = round(time_to_dg/30.4)) %>%
  group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo1) %>%
  slice(1) %>%
  ungroup()

df_3 = df_healthy %>%
  group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
  slice(2) %>%
  ungroup()
df_1 = df_healthy %>%
  dplyr::filter(!henkilotunnus %in% df_3$henkilotunnus) %>%
  group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
  slice(1) %>%
  ungroup() %>%
  dplyr::filter(min1 == 0)
df_healthy3 = rbind(df_3, df_1) %>%
  group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
  slice(1) %>%
  ungroup()
rm(df_1, df_3)

df_healthy1 = df_healthy3 %>%
  dplyr::arrange(desc(min1)) %>%
  group_by(henkilotunnus, time_to_dg_mo) %>%
  slice(1) %>%
  dplyr::mutate(tutkimus_lyhenne_yksikko = "Any") %>%
  ungroup()

df_healthy2 = df_healthy3 %>%
  dplyr::bind_rows(df_healthy1) %>%
  dplyr::mutate(disease = "healthy", state = 0, age = as.numeric(age)) %>%
  ungroup() %>%
  group_by(tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
  mutate(n = n(),
         prop = round(100 * sum(min1) / n, 1)) %>%
  ungroup()


# Load data
df2_1 = data.frame()
for (i in c("de novo AML", "MDS", "MF")) {
  
  print(i)
  
  # Load disease data
  df = readRDS(paste0(export, "/", i, "_lab_demo.rds")) %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% labtests) %>%
    dplyr::filter(time_to_dg >= (-5 * 365.25) & time_to_dg <= (-30)) %>%
    dplyr::mutate(time_to_dg_mo1 = round(time_to_dg/30.4)) %>%
    dplyr::mutate(min1 = ifelse(((tutkimus_lyhenne_yksikko == "B-Hb (g/l)" & sukupuoli_selite == "Nainen" & tulos < 120) |
                                   (tutkimus_lyhenne_yksikko == "B-Hb (g/l)" & sukupuoli_selite == "Mies" & tulos < 130) |
                                   (tutkimus_lyhenne_yksikko == "B-PLT (E9/l)" & tulos < 150) |
                                   (tutkimus_lyhenne_yksikko == "B-Neut (E9/l)" & tulos < 1.8)), 1, 0)) %>%
    # Select for each patients the lowest value per test and per year
    dplyr::arrange(desc(min1)) %>%
    # Max 1 sample per month
    dplyr::group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo1) %>%
    slice(1) %>%
    ungroup()
  
  df_3 = df %>%
    group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
    slice(2) %>%
    ungroup()
  df_1 = df %>%
    dplyr::filter(!henkilotunnus %in% df_3$henkilotunnus) %>%
    group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::filter(min1 == 0)
  df = rbind(df_3, df_1) %>%
    group_by(henkilotunnus, tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
    slice(1) %>%
    ungroup()
  rm(df_1, df_3)
  
  # Select for each patients the lowest value per year
  df1 = df %>%
    dplyr::arrange(desc(min1)) %>%
    group_by(henkilotunnus, time_to_dg_mo) %>%
    slice(1) %>%
    dplyr::mutate(tutkimus_lyhenne_yksikko = "Any") %>%
    ungroup()
  
  # Combine the two previous tables
  df2 = df %>%
    ungroup() %>%
    group_by(tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
    mutate(n = n(),
           n_max = mean(n)) %>%
    ungroup() %>%
    dplyr::bind_rows(df1) %>%
    group_by(time_to_dg_mo) %>%
    dplyr::mutate(
      n_min = mean(n_max, na.rm = TRUE),
      n = ifelse(tutkimus_lyhenne_yksikko == "Any", n_min, n)) %>%
    ungroup() %>%
    group_by(tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
    dplyr::mutate(prop = round(100*sum(min1, na.rm = TRUE) / n, 1),
                  state = 1,
                  age = as.numeric(age)) %>%
    ungroup()
  
  # Join
  df2_1 = rbind(df2_1, df2)
  
}

# Combine with similar tables for healthy
df3 = df2_1 %>%
  dplyr::bind_rows(df_healthy2)

df_healthy_sr = df3 %>%
  group_by(disease, tutkimus_lyhenne_yksikko) %>%
  summarise(q50 = mean(prop, na.rm=TRUE)) %>%
  ungroup(); df_healthy_sr

df_healthy_sr2 = df3 %>%
  dplyr::filter(disease == "healthy") %>%
  group_by(tutkimus_lyhenne_yksikko) %>%
  summarise(q50 = mean(prop, na.rm=TRUE),
            q50_low = q50-sd(prop, na.rm=TRUE),
            q50_high = q50+sd(prop, na.rm=TRUE)) %>%
  ungroup() %>%
  dplyr::select(-q50)
df_healthy_sr = df_healthy_sr %>%
  dplyr::left_join(df_healthy_sr2); df_healthy_sr

df2_1 = df2_1 %>%
  dplyr::left_join(df_healthy_sr)

df2_1 = df2_1 %>%
  dplyr::mutate(disease = ifelse(disease == "MF", "Primary MF", disease))


# Plot
for (i in unique(df2_1$tutkimus_lyhenne_yksikko)) {
  g = ggplot(df2_1 %>%
               dplyr::filter(tutkimus_lyhenne_yksikko == i), aes(x = time_to_dg_mo, y = prop)) +
    geom_rect(data = df2_1 %>%
                dplyr::filter(tutkimus_lyhenne_yksikko == i) %>%
                distinct(tutkimus_lyhenne, q50_low, q50_high, time_to_dg_mo, prop, .keep_all = FALSE),
              aes(ymin = q50_low, ymax = q50_high), xmin = -6, xmax = 0, fill="grey85") +
    geom_point(aes(color = disease), size = 3) +
    geom_line(aes(color = disease), size = 2) +
    labs(y="Proportion of patients with cytopenia (%)", x="Time to diagnosis (years)",
         title = gsub("Any", "Any cytopenia", i)) +
    scale_x_continuous(breaks = seq(from = -5, to = -1, by = 1), labels = seq(from = -5, to = -1, by = 1)) +
    ylim(0, 100) +
    theme_bw() +
    # scale_color_brewer(palette = "Set2") +
    scale_color_manual(values = c("#BF9F45","#348ABD", "#2B6E2A")) +
    guides(color=guide_legend(title="Disease", nrow = 1, override.aes = list(size = -1))) +
    theme(axis.line = element_line(colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position="bottom",
          legend.direction="horizontal",
          legend.title=element_text(colour="black", size=12),
          legend.title.align=0.5,
          legend.margin=margin(grid::unit(0, "cm")),
          legend.text=element_text(colour="black", size=12),
          legend.key.height=grid::unit(0.5, "cm"),
          legend.key.width=grid::unit(0.5, "cm"),
          plot.title.position = "panel",
          axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, colour = "black"),
          plot.title=element_text(colour="black", hjust = 0.0,
                                  # hjust=ifelse(i %in% c("secondary AML", "MF", "de novo AML"), 0.3, 0.4),
                                  vjust=0, size=14, face="bold")); g
  ggsave(plot = g,
         filename = paste0(results, "/cytopenia_proportion_", janitor::make_clean_names(i), ".png"),
         height = 5, width = 5, units = "in", dpi = 300)
}
