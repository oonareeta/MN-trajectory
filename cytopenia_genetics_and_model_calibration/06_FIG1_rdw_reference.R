# Heatmap

# Libraries
source("mounts/research/src/Rfunctions/library.R")
library(broom)
library(ggarchery)


# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Loop
for (i in c("MDS", "MF", "de novo AML")) {
  
  print(i)
  
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  # Remove f_p_kol_hdl, f_p_trfesat and f_p_fe_umol
  df1 = df1 %>%
    dplyr::filter(tutkimus_lyhenne_yksikko == "E-RDW (%)")
  
  # Categorical
  df1 = df1 %>%
    dplyr::mutate(tulos_cat = ifelse(sukupuoli_selite=="Mies" & tulos>=15 | sukupuoli_selite=="Nainen" & tulos>=14, 1, 0))
  
  # Summarise
  df2_low = df1 %>%
    arrange(tulos_cat) %>%
    group_by(henkilotunnus, time_to_dg_mo) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    group_by(time_to_dg_mo) %>%
    dplyr::summarise(
      tulos_low = 100*sum(tulos_cat, na.rm = TRUE) / n()
    ) %>%
    ungroup()
  df2_high = df1 %>%
    arrange(desc(tulos_cat)) %>%
    group_by(henkilotunnus, time_to_dg_mo) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    group_by(time_to_dg_mo) %>%
    dplyr::summarise(
      tulos_high = 100*sum(tulos_cat, na.rm = TRUE) / n()
    ) %>%
    ungroup()
  ## Combine
  df2 = full_join(df2_low, df2_high)
  ## Median
  df2$tulos = (df2$tulos_low+df2$tulos_high)/2

  print(paste0("The median of RDW is ",
               df2 %>%
                 filter(time_to_dg_mo > -5.1) %>%
                 group_by(time_to_dg_mo) %>%
                 summarise(median = round(median(tulos, na.rm=TRUE), 2))))

  # Colors
  if (i == "MDS") {
    cols1 = "#348ABD"
  } else if (i == "de novo AML") {
    cols1 = "#BF9F45"
  } else if (i == "MF") {
    cols1 = "#2B6E2A"
  }
  
  # Plot each regression with original values
  print("original")
  g = ggplot(df2[df2$time_to_dg_mo > -5.1 & df2$time_to_dg_mo < 0,], aes(x = time_to_dg_mo, y = tulos)) +
    geom_ribbon(aes(ymin = tulos_low, ymax = tulos_high), fill = cols1, alpha = 0.3) +
    geom_line() +
    geom_point() +
    labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
         y="E-RDW (%)", x=paste0("Time to diagnosis (years)")) +
    scale_x_continuous(breaks = seq(from = -5, to = -1, by = 1), labels = seq(from = -5, to = -1, by = 1)) +
    ylim(40, 100) +
    theme_bw() +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, colour = "black"),
          plot.title=element_text(size=14, face="bold", colour = "black"),
          axis.line = element_line(colour = "black"),
          title = element_text(size=14, face="bold", colour = "black"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank(),
          legend.position = "none")
  g
  # Export
  ggsave(plot = g,
         filename = paste0(results, "/RDW_reference_by_year_", i, ".png"),
         height = 5, width = 5, units = "in", dpi = 300)
}

