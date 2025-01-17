# Heatmap

# Author: Oscar Brück

# Libraries
source("./library.R")
library(broom)
library(ggarchery)


# General parameters
source("./parameters")


############################################################


# Data
## Read normalized healthy data
healthy = readRDS(paste0(import, "/lab_healthy_median.rds")) %>%
  dplyr::select(-n) %>%
  dplyr::rename(tulos_healthy = tulos)

# Loop
for (i in c("AML", "MDS", "MF", "de novo AML", "secondary AML")) {
  
  df = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  # Impact of age and gender
  df1 = df %>%
    dplyr::left_join(healthy)
  df1 = df1 %>%
    dplyr::mutate(tulos_norm = tulos - tulos_healthy)
  
  
  #  Collect significant lab and demographics data for healthy
  ## Adjust p
  ### Read linear regression tables - 5 y
  labtests5 = readxl::read_xlsx(paste0(results, "/", i, "/linear_regression_results_5y.xlsx"))
  labtests5$padj = p.adjust(labtests5$P_Value, method = "BH")
  labtests5$timetodg = -5
  ### Read linear regression tables - 4 y pre-disease
  labtests4 = readxl::read_xlsx(paste0(results, "/", i, "/linear_regression_results_4y.xlsx"))
  labtests4$padj = p.adjust(labtests4$P_Value, method = "BH")
  labtests4$timetodg = -4
  ### Read linear regression tables - 3 y pre-disease
  labtests3 = readxl::read_xlsx(paste0(results, "/", i, "/linear_regression_results_3y.xlsx"))
  labtests3$padj = p.adjust(labtests3$P_Value, method = "BH")
  labtests3$timetodg = -3
  ### Read linear regression tables - 2 y pre-disease
  labtests2 = readxl::read_xlsx(paste0(results, "/", i, "/linear_regression_results_2y.xlsx"))
  labtests2$padj = p.adjust(labtests2$P_Value, method = "BH")
  labtests2$timetodg = -2
  ### Read linear regression tables - 1 y pre-disease
  labtests1 = readxl::read_xlsx(paste0(results, "/", i, "/linear_regression_results_1y.xlsx"))
  labtests1$padj = p.adjust(labtests1$P_Value, method = "BH")
  labtests1$timetodg = -1
  
  # Join
  labtests = rbind(rbind(rbind(rbind(labtests5, labtests4), labtests3), labtests2), labtests1)
  
  labtests_2 = labtests2 %>%
    full_join(labtests3) %>%
    full_join(labtests4) %>%
    full_join(labtests5) %>%
    dplyr::arrange(padj) %>%
    dplyr::group_by(Variable) %>%
    dplyr::slice(1) %>%
    dplyr::ungroup() %>%
    distinct(Variable, .keep_all = TRUE) %>%
    dplyr::filter(!str_detect(Variable, "(E-Retic|E-RDW-SD|L-Myelos|L-Band|L-RBC|L-Metam|L-Segment|L-Erytrobl|AtypLy|L-Blast)"))
  
  labtests6 = labtests %>%
    dplyr::select(tutkimus_lyhenne_yksikko = Variable, time_to_dg_mo = timetodg, padj1_mo=padj)
  
  # Process data
  df2 = df1 %>%
    dplyr::filter(time_to_dg <(-30)) %>%
    dplyr::select(henkilotunnus, time_to_dg_mo, tutkimus_lyhenne_yksikko, tulos_norm) %>%
    dplyr::filter(tutkimus_lyhenne_yksikko %in% labtests_2$Variable) %>%
    group_by(time_to_dg_mo, tutkimus_lyhenne_yksikko) %>%
    summarise(tulos_norm = mean(tulos_norm, na.rm=TRUE),
              n = n()) %>%
    ungroup() %>%
    group_by(tutkimus_lyhenne_yksikko) %>%
    mutate(tulos_norm_mean = mean(tulos_norm, na.rm=TRUE),
           tulos_norm1 = round(tulos_norm-abs(tulos_norm_mean), 1),
           tulos_norm_scale = scales::rescale(tulos_norm1, to = c(0,1))) %>%
    ungroup() %>%
    dplyr::left_join(labtests_2 %>% dplyr::select(tutkimus_lyhenne_yksikko=Variable, padj, Coefficient)) %>%
    dplyr::left_join(labtests6) %>%
    dplyr::mutate(padj1 = ifelse(padj>=0.05, 0,
                                 ifelse(Coefficient<0, -1/padj, 1/padj)),
                  tutkimus_lyhenne_yksikko = gsub("^(L|B|fP)-", "", tutkimus_lyhenne_yksikko),
                  tutkimus_lyhenne_yksikko = ifelse(tutkimus_lyhenne_yksikko=="La (mm/h)", "ESR (mm/h)", tutkimus_lyhenne_yksikko),
                  size1 = ifelse(is.na(padj1_mo) | padj1_mo >= 0.05, 0, 1),
                  size = ifelse(is.na(padj1_mo) | padj1_mo >= 0.05, 0,
                                ifelse(padj1_mo<0.001, 3,
                                       ifelse(padj1_mo<0.01, 2, 1))))
  
  tt_5_2 = df2 %>%
    dplyr::filter(time_to_dg_mo > -3.1 & time_to_dg_mo < (-1)) %>%
    group_by(tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
    summarise(tulos_norm_scale = mean(tulos_norm_scale, na.rm=TRUE)) %>%
    ungroup() %>%
    group_by(tutkimus_lyhenne_yksikko) %>%
    summarise(tulos_norm_scale_pre = mean(tulos_norm_scale, na.rm=TRUE)) %>%
    ungroup()
  tt_1 = df2 %>%
    dplyr::filter(time_to_dg_mo == -1) %>%
    group_by(tutkimus_lyhenne_yksikko, time_to_dg_mo) %>%
    summarise(tulos_norm_scale = mean(tulos_norm_scale, na.rm=TRUE)) %>%
    ungroup() %>%
    group_by(tutkimus_lyhenne_yksikko) %>%
    summarise(tulos_norm_scale_post = mean(tulos_norm_scale, na.rm=TRUE)) %>%
    ungroup()
  tt_5_2 = full_join(tt_5_2, tt_1) %>%
    dplyr::mutate(padj2 = tulos_norm_scale_post-0.5) %>%
    dplyr::select(tutkimus_lyhenne_yksikko, padj2)
  df2 = df2 %>%
    dplyr::left_join(tt_5_2) %>%
    dplyr::mutate(padj1 = ifelse((padj1<0 & padj2>0) | (padj1>0 & padj2<0), -padj1, padj1))
  
  # Plot
  y1 = df2 %>%
    dplyr::filter(time_to_dg_mo == -5) %>%
    distinct(tutkimus_lyhenne_yksikko, .keep_all = TRUE) %>%
    dplyr::filter(!is.na(padj1) & padj1 < 0) %>%
    nrow(); y1
  y2 = df2 %>%
    dplyr::filter(time_to_dg_mo == -5) %>%
    dplyr::filter(!is.na(padj1) & padj1 > 0) %>%
    nrow(); y2
  y3 = df2 %>%
    dplyr::filter(time_to_dg_mo == -5) %>%
    dplyr::filter(!is.na(padj1) & padj1 == 0) %>%
    nrow(); y3
  
  g = ggplot(df2 %>%
               dplyr::filter(time_to_dg_mo > -5.1),
             aes(x=time_to_dg_mo, y=reorder(tutkimus_lyhenne_yksikko, padj1), fill=tulos_norm_scale)) +
    geom_tile(colour="white", size=0.5) +
    geom_arrowsegment(aes(x = 0, y = y1, xend = 0, yend = 1),
                      fill = "blue", color = "blue", size = 0.5, 
                      arrows = arrow(type = "closed", length = unit(0.15, "inches"))) +
    geom_arrowsegment(aes(x = 0, y = (y1+1+y3), xend = 0, yend = (y1+y2+y3)), 
                      fill = "red", color = "red", size = 0.5, 
                      arrows = arrow(type = 'closed', length = unit(0.15, "inches"))) +
    annotate("text", x=0.5, y=y1/2, label=paste0("AdjP", sprintf('\u2193')), fontface="bold", angle=90) + 
    annotate("text", x=0.5, y=((y1+y3)+(y2/2)), label=paste0("AdjP", sprintf('\u2193')), fontface="bold", angle=90) + 
    geom_point(aes(size = size1), color="black") +
    guides(fill=guide_legend(title="Normalized\nvalues", override.aes = list(size = -1)),
           size=guide_legend(title="AdjP")) +
    labs(x="Time to diagnosis (years)", y="", title = ifelse(i=="MF", "Primary MF",
                                                             ifelse(i=="de novo AML", "De novo AML",
                                                                    ifelse(i=="secondary AML", "Secondary AML", i)))) +
    scale_fill_distiller(palette = "RdBu") +
    scale_size_continuous(range = c(0, 2),
                          breaks = c(1),
                          limits = c(0.1, 2),
                          labels = c("p<0.05")) +
    scale_x_continuous(breaks = seq(from = -5, to = -1, by = 1),
                       labels = seq(from = -5, to = -1, by = 1)) +#labels = seq(from = -60, to = 0, by = 12)) +
    coord_fixed() +
    theme_bw() +
    theme(legend.position="right",
          legend.direction="vertical",
          legend.title=element_text(colour="black"),
          legend.title.align=0.5,
          legend.margin=margin(grid::unit(0, "cm")),
          legend.text=element_text(colour="black", size=7, face="bold"),
          legend.key.height=grid::unit(0.5, "cm"),
          legend.key.width=grid::unit(0.5, "cm"),
          plot.title.position = "panel",
          axis.text.x=element_text(size=10, colour="black"),
          axis.text.y=element_text(colour="black"),
          #axis.title.x =element_text(colour="black", hjust=1),
          axis.ticks=element_line(size=0.5),
          # plot.background=element_blank(),
          plot.margin=margin(0.7, 0.4, 0.1, 0.2, "cm"),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          legend.box.margin=margin(-5,-5,-5,-5),
          plot.title=element_text(colour="black",
                                  hjust=ifelse(i %in% c("secondary AML", "MF", "de novo AML"), 0.3, 0.4),
                                  vjust=-1, size=12, face="bold")
    ); g
  ## Export
  ggsave(plot = g,
         filename = paste0(results, "/", i, "/heatmap.png"),
         height = 6, width = 6, units = "in", dpi = 300)
}
