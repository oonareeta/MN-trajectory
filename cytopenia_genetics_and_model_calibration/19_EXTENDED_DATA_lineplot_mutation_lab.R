# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")


# Load test data
disease = c("de_novo_AML", "MDS", "primary_MF")

# Prediction probability by mutations and karyotype alterations
# Loop over diseases
for (i in disease) {
  
  # Load test data
  df = fread(paste0("mounts/research/husdatalake/disease/processed_data/Preleukemia/lagged_data_", gsub("_", " ", gsub("primary_", "", i)), "3.csv"),
             select = c("henkilotunnus", "time_to_dg", "disease", "b_baso_e9_l_tulos_norm", "b_eos_e9_l_tulos_norm", "b_erblast_e9_l_tulos_norm", "b_eryt_e12_l_tulos_norm", "b_hb_g_l_tulos_norm",
                        "b_hkr_percent_tulos_norm", "b_ly_e9_l_tulos_norm", "b_monos_e9_l_tulos_norm", "b_neut_e9_l_tulos_norm", "b_plt_e9_l_tulos_norm",
                        "b_wbc_e9_l_tulos_norm", "e_mch_pg_tulos_norm", "e_mchc_g_l_tulos_norm", "e_mcv_fl_tulos_norm", "e_rdw_percent_tulos_norm",
                        "e_retic_percent_tulos_norm", "l_baso_percent_tulos_norm", "l_eos_percent_tulos_norm", "l_ly_percent_tulos_norm",
                        "l_mono_percent_tulos_norm", "l_neut_percent_tulos_norm"))
  
  # Folder
  if (i == "primary_MF") {
    j = "MF"
  } else if (i == "de_novo_AML") {
    j = "de novo AML"
  } else {
    j = i
  }
  
  dir.create(paste0(results, "/", j, "/Alteration_lab_lineplot"))
  
  if (i == "any_MN") {
    # Load cytogenetics and mutation data
    df_c = fread(paste0(export, "/chromosomes_lab_MDS.csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct() %>%
      full_join(fread(paste0(export, "/chromosomes_lab_de novo AML.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct()) %>%
      full_join(fread(paste0(export, "/chromosomes_lab_MF.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct())
    df_m = fread(paste0(export, "/mutations_lab_MDS.csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct() %>%
      full_join(fread(paste0(export, "/mutations_lab_de novo AML.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct()) %>%
      full_join(fread(paste0(export, "/mutations_lab_MF.csv")) %>%
                  dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
                  distinct())
  } else {
    # Load cytogenetics and mutation data
    df_c = fread(paste0(export, "/chromosomes_lab_", j, ".csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct()
    df_m = fread(paste0(export, "/mutations_lab_", j, ".csv")) %>%
      dplyr::select(-c(time_to_dg, tutkimus_lyhenne_yksikko, tulos, tulos_healthy, tulos_norm, sukupuoli_selite, age)) %>%
      distinct()
  }
  
  # Remove rows where all variables from 2 to ncol are ""
  df_c <- df_c[!apply(df_c[, 2:ncol(df_c)], 1, function(row) all(row == "")), ]
  df_m <- df_m[!apply(df_m[, 2:ncol(df_m)], 1, function(row) all(row == "")), ]
  
  
  # Names of alterations
  chrm_pos_counts = names(df_c)[2:ncol(df_c)]
  mut_pos_counts = names(df_m)[2:ncol(df_m)]
  alterations = c(chrm_pos_counts, mut_pos_counts)
  labs = c("b_baso_e9_l_tulos_norm", "b_eos_e9_l_tulos_norm", "b_erblast_e9_l_tulos_norm", "b_eryt_e12_l_tulos_norm", "b_hb_g_l_tulos_norm",
           "b_hkr_percent_tulos_norm", "b_ly_e9_l_tulos_norm", "b_monos_e9_l_tulos_norm", "b_neut_e9_l_tulos_norm", "b_plt_e9_l_tulos_norm",
           "b_wbc_e9_l_tulos_norm", "e_mch_pg_tulos_norm", "e_mchc_g_l_tulos_norm", "e_mcv_fl_tulos_norm", "e_rdw_percent_tulos_norm",
           "e_retic_percent_tulos_norm", "l_baso_percent_tulos_norm", "l_eos_percent_tulos_norm", "l_ly_percent_tulos_norm",
           "l_mono_percent_tulos_norm", "l_neut_percent_tulos_norm")
  names1 = readxl::read_xlsx("mounts/research/husdatalake/disease/processed_data/Preleukemia/variable_names.xlsx")
  
  
  # Join mutations and chromosomes
  df_c = df_c %>%
    dplyr::inner_join(df)
  df_m = df_m %>%
    dplyr::inner_join(df)
  df1 = full_join(df_c, df_m) %>%
    mutate(time_to_dg = -time_to_dg / 365.24)
    # mutate(time_to_dg = -round(time_to_dg / 365.24))
  
  for (alteration1 in alterations) {
    
    print(alteration1)
    
    for (lab1 in labs) {
      
      df2 = df1
      df2$x1 = df2[[alteration1]]
      df2$y1 = df2[[lab1]]
      
      y_lab = paste0(names1[names1$variable==lab1,]$output, ", normalized", sep="")
      
      # Legend title
      legend1 =
        gsub("_", ",",
             gsub("trisomy", "Trisomy ",
                  gsub("monosomy", "Monosomy ",
                       gsub("monosomy_y", "Monosomy Y",
                            # gsub("(monosomy|monosomy_)", "-",
                            gsub("_kar", " kar",
                                 
                                 gsub("karyotype_HR", "High-risk karyotype",
                                      gsub("translocation", "t(",
                                           gsub("t_", "t(",
                                                gsub("inv_", "inv(",
                                                     gsub("none", "Normal karyotype",
                                                          gsub("del_", "del(", alteration1)))))))))))
      legend1 = ifelse(str_detect(legend1, "(^(t\\(|inv|del))"), paste0(legend1, ")"), legend1)
      legend1 = gsub("mono", "Mono",
                     gsub("tri", "Tri",
                          gsub("com", "Com",
                               gsub("t\\(3\\)", "t(3,)",
                                    gsub("t\\(11\\)", "t(11,)", legend1)))))
      legend1 = ifelse(legend1 == "None", "No mutation", gsub("FLT3,", "FLT3-", legend1))
      
      
      # Remove
      df2 = df2 %>%
        dplyr::filter(!x1 == "")
      
      # Plot
      # b_neut_e9_l_tulos_norm
      g = ggplot(df2, aes(x = time_to_dg, y = y1, color = x1)) +
        geom_smooth(method = "loess", se = TRUE, fill = "grey80") + #method = "loess",
        labs(title = paste0(toupper(substr(gsub("_", " ", i), 1, 1)), substr(gsub("_", " ", i), 2, nchar(gsub("_", " ", i))), sep=""),
             y=y_lab,
             x="Time to diagnosis (years)") +
        scale_color_brewer(name = legend1, palette = "Set1", direction = -1) +
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
              legend.position = "bottom"); g
      
      # Export
      ggsave(plot = g,
             filename = paste0(results, "/", j, "/Alteration_lab_lineplot/Lineplot_", alteration1, "_", lab1, ".png"),
             height = 4, width = 4, units = "in", dpi = 300)
      
    }
    
  }
  
}
