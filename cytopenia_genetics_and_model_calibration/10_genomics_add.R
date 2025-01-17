# Link cytogenetics and genomics data to prediction accuracy
## The idea is to estimate how much any given mutation affects pre-leukemia lab values

# Author: Oscar Brück

# Libraries
source("./library.R")
library(cowplot)


# General parameters
source("./parameters")


######## Prepare data ########


# Loop over diseases
diseases = c("MDS", "MF", "de novo AML") #c("AML", "MDS", "MF", "de novo AML", "secondary AML")
for (i in diseases) {
  
  dir.create(paste0(results, "/", i, "/KM"))
  dir.create(paste0(results, "/", i, "/Scatterplots"))
  
  # Collect data for disease
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  # Read lab data
  df = fread(paste0(export, "/lab_data_for_modelling_", i, ".csv")) %>%
    dplyr::filter(henkilotunnus %in% unique(df1$henkilotunnus)) %>%
    dplyr::filter(time_to_dg<-30 & time_to_dg>-730)
  
  # La to ESR
  df = df %>%
    dplyr::mutate(tutkimus_lyhenne_yksikko = ifelse(tutkimus_lyhenne_yksikko=="B-La (mm/h)", "B-ESR (mm/h)", tutkimus_lyhenne_yksikko))
  
  # Summarise median by patient and laboratory test
  df = df %>%
    dplyr::group_by(henkilotunnus, tutkimus_lyhenne_yksikko) %>%
    mutate(tulos_norm = mean(tulos_norm, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup()
  
  
  # Molecular genetics
  if (i %in% c("AML", "de novo AML", "secondary AML")) {
    mut = readRDS("mounts/research/husdatalake/disease/processed_data/AML/ELN.rds") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
    chrms = mut %>%
      dplyr::select(which(names(mut) == "henkilotunnus") : ((which(names(mut) == "ABL1")-1)), translocation3, translocation11, complex_karyotype, monosomal_karyotype, karyotype_HR) %>%
      dplyr::mutate_at(.vars = c("translocation3", "translocation11", "complex_karyotype", "monosomal_karyotype", "karyotype_HR"), function(x) ifelse(is.na(x), NA,
                                                                                                                                                      ifelse(x==1, "POS", "NEG")))
    pat = fread("mounts/research/husdatalake/disease/processed_data/AML/AML_pathology.csv") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus) %>%
      dplyr::mutate(mgg_value_percentage = ifelse(mgg_value_percentage > 100, mgg_value_percentage/100, mgg_value_percentage)) %>%
      dplyr::select(henkilotunnus, flow_blast_percentage, mgg_value_percentage) %>%
      dplyr::mutate(bm_blast = pmax(flow_blast_percentage, mgg_value_percentage, na.rm=TRUE)) %>%
      dplyr::filter(bm_blast > 9.9)
    
    tr3 = fread("mounts/research/husdatalake/disease/processed_data/AML/AML_survival_table.csv")
    tr4 = tr3 %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
  } else if (i %in% c("MDS")) {
    mut = readRDS("mounts/research/husdatalake/disease/processed_data/MDS/IPSSR.rds") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
    chrms = mut %>%
      dplyr::select(which(names(mut) == "henkilotunnus") : ((which(names(mut) == "ABL1")-1)), translocation3, translocation11, complex_karyotype, monosomal_karyotype, karyotype_HR) %>%
      dplyr::mutate_at(.vars = c("translocation3", "translocation11", "complex_karyotype", "monosomal_karyotype", "karyotype_HR"), function(x) ifelse(is.na(x), NA,
                                                                                                                                                      ifelse(x==1, "POS", "NEG")))
    pat = fread("mounts/research/husdatalake/disease/processed_data/MDS/etiology.csv") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus) %>%
      dplyr::select(henkilotunnus, bm_blast = bm_blast_mgg)
    
    tr3 = fread("mounts/research/husdatalake/disease/processed_data/MDS/demographics.csv") %>%
      dplyr::mutate(OS_event = ifelse(!is.na(kuolinaika_pvm), 1, 0),
                    OS_time = ifelse(OS_event == 1, (kuolinaika_pvm - dg_date_combined) / 365.24, as.numeric(as.Date(gsub("_", "/", read.table("mounts/research/foldername")$V1)) - as.Date(dg_date_combined)) / 365.24)) %>%
      dplyr::select(henkilotunnus, OS_event, OS_time)
    
    tr4 = tr3 %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
  } else {
    mut = readRDS("mounts/research/husdatalake/disease/processed_data/MF/ELN.rds") %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
    
    chrms = mut %>%
      dplyr::select(which(names(mut) == "henkilotunnus") : ((which(names(mut) == "ABL1")-1)), complex_karyotype, monosomal_karyotype) %>%
      dplyr::mutate_at(.vars = c("complex_karyotype", "monosomal_karyotype"), function(x) ifelse(is.na(x), NA,
                                                                                                 ifelse(x==1, "POS", "NEG")))
    
    tr3 = readRDS("mounts/research/husdatalake/disease/processed_data/MF/demographics.rds") %>%
      dplyr::mutate(OS_event = ifelse(!is.na(kuolinaika_pvm), 1, 0),
                    OS_time = ifelse(OS_event == 1, (as.Date(kuolinaika_pvm) - as.Date(dg_date_combined)) / 365.24, as.numeric(as.Date(gsub("_", "/", read.table("mounts/research/foldername")$V1)) - as.Date(dg_date_combined)) / 365.24)) %>%
      dplyr::select(henkilotunnus, OS_event, OS_time)
    
    tr4 = tr3 %>%
      dplyr::filter(henkilotunnus %in% df$henkilotunnus)
  }
  
  
  # Edit none
  mut = mut %>%
    dplyr::mutate(None = ifelse(is.na(None) & !is.na(ASXL1), "NEG", None))
  
  
  # Mutation data
  muts = mut %>%
    dplyr::select(henkilotunnus, which(names(mut) == "ABL1") : which(names(mut) == "ZRSR2"))
  
  
  ## Keep columns with ≥10 POS
  # Function to count POS per column
  count_pos_per_column <- function(df) {
    sapply(df, function(x) sum(x == "POS", na.rm=TRUE))
  }
  
  # Count NAs in the sample data frame
  if (i == "MF") {
    n1 = 4
  } else {
    n1 = 10
  }
  
  ## Chromosomes
  chrm_pos_counts = count_pos_per_column(chrms)
  chrm_pos_counts = chrm_pos_counts[chrm_pos_counts>=n1]
  chrm_pos_counts = chrm_pos_counts[names(chrm_pos_counts)[grep(pattern = "der_", x = names(chrm_pos_counts), invert = TRUE)]]
  ## Mutations
  mut_pos_counts = count_pos_per_column(muts)
  mut_pos_counts = mut_pos_counts[mut_pos_counts>=n1]
  
  # Combine
  df_c = df %>%
    dplyr::inner_join(chrms %>% dplyr::select(henkilotunnus, all_of(names(chrm_pos_counts))))
  df_c = df_c %>%
    dplyr::select(-contains("t_11q_3"), -starts_with("der")) %>%
    dplyr::left_join(tr3 %>%
                       dplyr::select(henkilotunnus, OS_time, OS_event)) %>%
    dplyr::left_join(pat %>%
                       dplyr::select(henkilotunnus, bm_blast))
  df_m = df %>%
    dplyr::inner_join(muts %>% dplyr::select(henkilotunnus, all_of(names(mut_pos_counts)))) %>%
    dplyr::left_join(tr3 %>%
                       dplyr::select(henkilotunnus, OS_time, OS_event)) %>%
    dplyr::left_join(pat %>%
                       dplyr::select(henkilotunnus, bm_blast))
  
  if (!i == "MF") {
    df_c = df_c %>%
      dplyr::left_join(pat %>%
                         dplyr::select(henkilotunnus, bm_blast))
    
    df_m = df_m %>%
      dplyr::left_join(pat %>%
                         dplyr::select(henkilotunnus, bm_blast))
  }
  
  if (i == "de novo AML") {
    df_c = df_c %>%
      dplyr::left_join(tr3 %>%
                         dplyr::select(henkilotunnus, HDC))
    
    df_m = df_m %>%
      dplyr::left_join(tr3 %>%
                         dplyr::select(henkilotunnus, HDC))
  }
  
  
  ############## WILCOXON TEST ##############
  
  
  for (j in unique(df_m$tutkimus_lyhenne_yksikko)) {
    for (k in (names(df_m)[((which(names(df_m) == "age")+1) : (which(names(df_m) == "OS_time")-1))])) {
      
      df_m1 = df_m %>%
        dplyr::filter(tutkimus_lyhenne_yksikko == j)
      df_m1$gene = df_m1[[k]]
      
      df_m1 = df_m1 %>%
        dplyr::filter(!is.na(gene)) %>%
        dplyr::mutate(gene = ifelse(gene == "POS", "MUT", ifelse(gene == "NEG", "WT", NA)),
                      gene = factor(gene, levels = c("WT", "MUT")))
      
      g = ggplot(data = df_m1, aes(x = gene, y = tulos_norm)) + #, group=time_to_dg_mo
        geom_jitter(width = 0.2, size = 0.5, color = "black") +
        geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
        labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
             y=j, x=k) +
        theme_bw() +
        stat_compare_means(method = "wilcox.test", label.x = 1.25) +#, aes(label = ..p.signif..)) +
        theme(axis.text.x = element_text(size=12, colour = "black"),
              axis.text.y = element_text(size=12, colour = "black"),
              axis.title=element_text(size=14, colour = "black"),
              plot.title=element_text(size=14, face="bold", colour = "black"),
              legend.position = "none",
              panel.border = element_blank(),
              axis.line = element_line(colour = "black"),
              panel.background = element_blank(),
              panel.grid.minor = element_blank()); g
      ggsave(plot = g,
             filename = paste0(results, "/", i, "/Scatterplots/", janitor::make_clean_names(j), "_", janitor::make_clean_names(k), "_norm.png"),
             height = 4, width = 4, units = "in", dpi = 300)
      
    }
  }

  
  # Only plot for TP53
  if (i %in% c("MDS", "de novo AML")) {
    tmp1 = df_m %>%
      dplyr::filter(!is.na(TP53)) %>%
      dplyr::mutate(TP53 = ifelse(TP53 == "POS", "MUT", ifelse(TP53 == "NEG", "WT", NA)),
                    TP53 = factor(TP53, levels = c("WT", "MUT"))) %>%
      dplyr::filter(tutkimus_lyhenne_yksikko == "E-RDW (%)")
    
    g = ggplot(data = tmp1, aes(x = TP53, y = tulos_norm)) + #, group=time_to_dg_mo
      geom_jitter(width = 0.2, size = 0.5, color = "black") +
      geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
      labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
           y="E-RDW (%)", x="TP53") +
      theme_bw() +
      stat_compare_means(method = "wilcox.test", label.x = 1.25) +#, aes(label = ..p.signif..)) +
      theme(axis.text.x = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title=element_text(size=14, colour = "black"),
            plot.title=element_text(size=14, face="bold", colour = "black"),
            legend.position = "none",
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.background = element_blank(),
            panel.grid.minor = element_blank()); g
    ggsave(plot = g,
           filename = paste0(results, "/", i, "/RDW_TP53_norm.png"),
           height = 4, width = 4, units = "in", dpi = 300)
    
    g = ggplot(data = tmp1, aes(x = TP53, y = tulos)) + #, group=time_to_dg_mo
      geom_jitter(width = 0.2, size = 0.5, color = "black") +
      geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
      labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
           y="E-RDW (%)", x="TP53") +
      theme_bw() +
      stat_compare_means(method = "wilcox.test", label.x = 1.25) +#, aes(label = ..p.signif..)) +
      theme(axis.text.x = element_text(size=12, colour = "black"),
            axis.text.y = element_text(size=12, colour = "black"),
            axis.title=element_text(size=14, colour = "black"),
            plot.title=element_text(size=14, face="bold", colour = "black"),
            legend.position = "none",
            panel.border = element_blank(),
            axis.line = element_line(colour = "black"),
            panel.background = element_blank(),
            panel.grid.minor = element_blank()); g
    ggsave(plot = g,
           filename = paste0(results, "/", i, "/RDW_TP53.png"),
           height = 4, width = 4, units = "in", dpi = 300)
  }
  
  
  ############## COX REGRESSION TESTS ##############
  
  
  # Loop karyotypes
  res <- data.frame()
  for (lab1 in 1:length(unique(df$tutkimus_lyhenne_yksikko))) {
    print(paste0(lab1, "/", length(unique(df$tutkimus_lyhenne_yksikko))))
    
    df_c1 = df_c %>%
      dplyr::filter(tutkimus_lyhenne_yksikko == unique(df$tutkimus_lyhenne_yksikko)[lab1])# %>%
    
    if (i == "de novo AML") {
      df_c1 = df_c %>%
        dplyr::filter(tutkimus_lyhenne_yksikko == unique(df$tutkimus_lyhenne_yksikko)[lab1]) %>%
        group_by(henkilotunnus, OS_time, OS_event, HDC) %>%
        summarise(tulos_norm = mean(tulos_norm, na.rm = TRUE),
                  tulos = mean(tulos, na.rm = TRUE)) %>%
        ungroup() %>% dplyr::filter(HDC == 1)
      
      print(paste0("There are ", length(unique(df_c$henkilotunnus)), " AML patients of which ", length(unique(df_c1$henkilotunnus)), " are treated intensively."))
    }
    
    df_c1 = df_c1 %>%
      dplyr::mutate(tulos_norm_scaled = scales::rescale(x = tulos_norm, to = c(0,1))) %>%
      dplyr::mutate(tulos_norm_cat = ntile(tulos_norm, 3))
    
    # Cox regression
    multiple_t_tests_p_value = summary(coxph(Surv(OS_time, OS_event)~tulos_norm_scaled, data = df_c1))
    pvalue <- signif(multiple_t_tests_p_value$coefficients[5], digits=5)
    padj <- pvalue
    beta <- signif(multiple_t_tests_p_value$coef[1], digits=5);#coeficient beta
    HR <- signif(multiple_t_tests_p_value$coef[2], digits=5);#exp(beta)
    HR.confint.lower <- signif(multiple_t_tests_p_value$conf.int[1,"lower .95"], 5)
    HR.confint.upper <- signif(multiple_t_tests_p_value$conf.int[1,"upper .95"], 5)
    res1 <- c(unique(df$tutkimus_lyhenne_yksikko)[lab1], i, beta, HR, HR.confint.lower, HR.confint.upper, pvalue, padj)
    rm(beta, HR, HR.confint.lower, HR.confint.upper, pvalue, padj)
    names(res1) <- c("covariate", "Disease", "beta", "HR", "CI95low", "CI95high", "pvalue", "padj")
    
    # Save data
    res <- rbind(res, res1)
    names(res)<-c("covariate", "Disease", "beta", "HR", "CI95low", "CI95high", "pvalue", "padj")
    
    # Adjust p value
    res$padj <- p.adjust(res$pvalue, method = "BH")
    
    # Plot
    fit = survfit(Surv(OS_time, OS_event)~tulos_norm_cat, data = df_c1)
    
    if (length(fit$n) > 1) {
      g = ggsurvplot(fit,
                     data = df_c1,
                     palette = "Set1",
                     size = 2,   #line thickness
                     ggtheme = theme_minimal(), #theme
                     font.main = c(15, "black"), #title font
                     font.x = c(15, "bold", "black"), #x-axis font
                     font.y = c(15, "bold", "black"), #y-axis font
                     font.tickslab = c(15, "bold", "black"), #axis numbering font
                     conf.int = FALSE, #confidence interval
                     pval = TRUE, #p-value
                     pval.size = 5, #p-value size
                     pval.coord = c(0, 0.1),
                     risk.table.pos = "out",
                     tables.y.text = FALSE,
                     tables.theme = theme_cleantable(),
                     break.x.by = 2,
                     risk.table.fontsize  = 4,
                     risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                     risk.table.col = "strata", #risk table color
                     risk.table.height = 0.25,
                     ylab = "OS event (%)",
                     xlab = "Time (years)",
                     censor = TRUE,
                     censor.shape = 108,
                     censor.size = 4,
                     font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
                     legend.title = unique(df$tutkimus_lyhenne_yksikko)[lab1],
                     legend.labs = c("Low", "Intermediate", "High"))
      
      ggsave(plot = g, filename = paste0(results, "/", i, "/KM/", janitor::make_clean_names(unique(df$tutkimus_lyhenne_yksikko)[lab1]), ".png"), width = 6, height = 6, bg = "white", units = 'in', dpi = 300)
    } 
  }
  
  # Export
  writexl::write_xlsx(res, paste0(results, "/", i, "/KM/Cox_regression_univariate.xlsx"))
  
}


############## VOLCANO PLOT ##############


res = data.frame()
for (i in diseases) {
  res = rbind(res, readxl::read_xlsx(paste0(results, "/", i, "/KM/Cox_regression_univariate.xlsx")))
}

# Prepare data for volcanoplot
res1 <- res %>% mutate(
  fold.change.log10 = log(as.numeric(HR), 10),
  pvalue = as.numeric(pvalue),
  sig = ifelse(pvalue>=0.05, "ns",
               ifelse(HR >= 1, "HR≥1", "HR<1"))
)

# Plot
g = ggplot(res1, aes(fold.change.log10, -log10(pvalue))) +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "black", alpha = 0.6) +
  geom_point(aes(fill=Disease, size=-log10(pvalue)), shape = 21, color="black") +
  labs(x="Log10 Hazard ratio", y="Log10 p-value") +
  guides(fill=guide_legend("Disease",
                           # title.position = "top",
                           title.hjust = 0.5,
                           override.aes = list(size = 5)),
         size = FALSE) +
  theme_bw() +
  xlim(-1, 1) +
  theme(legend.direction = 'horizontal', 
        legend.key = element_rect(size = 5),
        legend.key.size = unit(1.5, 'lines'),
        axis.text.x = element_text(size=12, colour = "black"),
        axis.text.y = element_text(size=12, colour = "black"),
        axis.title=element_text(size=14, colour = "black"),
        plot.title=element_text(size=14, face="bold", colour = "black"),
        legend.title=element_text(size=12, face="bold", colour = "black"),
        legend.text=element_text(size=12, colour = "black"),
        panel.grid.minor = element_blank(),
        legend.position = "bottom") +
  scale_fill_manual(values=c("#BF9F45", "#348ABD", "#2B6E2A")) +
  geom_text_repel(data=res1[res1$pvalue<0.05,], size=4, aes(label=covariate), box.padding = 0.5);g
ggsave(plot = g, filename = paste0(results, "/Volcano_OS_labs.png"), width = 6, height = 6, dpi = 300, units = "in")


############## RDW AT DIAGNOSIS ##############


for (i in diseases) {
  
  # Read RDW data
  lab_dg = readRDS(paste0("mounts/research/husdatalake/disease/processed_data/", gsub("de novo ", "", i), "/lab_dg.rds")) %>%
    dplyr::filter(str_detect(tutkimus_lyhenne, "E -RDW")) %>%
    dplyr::select(henkilotunnus, RDW=tulos_mean)
  
  # Read survival data
  if (i == "de novo AML") {
    
    tr3 = fread("mounts/research/husdatalake/disease/processed_data/AML/AML_survival_table.csv")
    
  } else if (i == "MDS") {
    
    tr3 = fread("mounts/research/husdatalake/disease/processed_data/MDS/demographics.csv") %>%
      dplyr::mutate(OS_event = ifelse(!is.na(kuolinaika_pvm), 1, 0),
                    OS_time = ifelse(OS_event == 1, (kuolinaika_pvm - dg_date_combined) / 365.24, as.numeric(as.Date(gsub("_", "/", read.table("mounts/research/foldername")$V1)) - as.Date(dg_date_combined)) / 365.24)) %>%
      dplyr::select(henkilotunnus, OS_event, OS_time)
    
  } else if (i == "MF") {
    
    tr3 = readRDS("mounts/research/husdatalake/disease/processed_data/MF/demographics.rds") %>%
      dplyr::mutate(OS_event = ifelse(!is.na(kuolinaika_pvm), 1, 0),
                    OS_time = ifelse(OS_event == 1, (as.Date(kuolinaika_pvm) - as.Date(dg_date_combined)) / 365.24, as.numeric(as.Date(gsub("_", "/", read.table("mounts/research/foldername")$V1)) - as.Date(dg_date_combined)) / 365.24)) %>%
      dplyr::select(henkilotunnus, OS_event, OS_time)
    
    # Remove secondary MF
    pmf = readRDS("mounts/research/husdatalake/disease/processed_data/MF/secondary_MF.rds") %>%
      dplyr::filter(disease == "Primary")
    tr3 = tr3 %>%
      dplyr::filter(henkilotunnus %in% pmf$henkilotunnus)
    
  }
  
  ## Join
  lab_dg1 = lab_dg %>%
    dplyr::inner_join(tr3)
  
  ## Categorize
  lab_dg1$RDW1 = ntile(lab_dg1$RDW, 3)
  
  lab_dg1 %>% group_by(RDW1) %>% summarise(RDW = mean(RDW), RDW_min = min(RDW))
  
  fit = survfit(Surv(OS_time, OS_event)~RDW1, data = lab_dg1)
  
  g = ggsurvplot(fit,
                 data = lab_dg1,
                 palette = "Set1",
                 size = 2,   #line thickness
                 ggtheme = theme_minimal(), #theme
                 font.main = c(15, "black"), #title font
                 font.x = c(15, "bold", "black"), #x-axis font
                 font.y = c(15, "bold", "black"), #y-axis font
                 font.tickslab = c(15, "bold", "black"), #axis numbering font
                 conf.int = FALSE, #confidence interval
                 pval = TRUE, #p-value
                 pval.size = 5, #p-value size
                 pval.coord = c(0, 0.1),
                 risk.table.pos = "out",
                 tables.y.text = FALSE,
                 tables.theme = theme_cleantable(),
                 break.x.by = 2,
                 risk.table.fontsize  = 4,
                 risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                 risk.table.col = "strata", #risk table color
                 risk.table.height = 0.25,
                 ylab = "OS event (%)",
                 xlab = "Time (years)"
                 censor = TRUE,
                 censor.shape = 108,
                 censor.size = 4,
                 font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
                 legend.title = "E-RDW (%)",
                 legend.labs = c("Low", "Intermediate", "High")); g
  ggsave(plot = g, filename = paste0(results, "/", i, "/KM/RDW_OS_all_patients.png"), width = 6, height = 6, bg = "white", units = 'in', dpi = 300)
  
  
  
  if (i %in% c("MDS", "de novo AML")) {
    
    # Plot
    if (i == "MDS") {
      eln = readRDS("mounts/research/husdatalake/disease/processed_data/MDS/IPSSR.rds") %>%
        dplyr::select(henkilotunnus, TP53)
    } else {
      eln = readRDS("mounts/research/husdatalake/disease/processed_data/AML/ELN.rds") %>%
        dplyr::select(henkilotunnus, TP53)
    }
    tr4 = tr3 %>%
      dplyr::inner_join(eln)
    
    lab_dg2 = full_join(tr4 %>% dplyr::select(henkilotunnus, OS_event, OS_time, TP53),
                        lab_dg1 %>% dplyr::select(henkilotunnus, OS_event, OS_time, RDW=RDW1))
    
    fit1 = survfit(Surv(OS_time, OS_event)~TP53, data = lab_dg2)
    fit2 = survfit(Surv(OS_time, OS_event)~RDW, data = lab_dg2)
    fit <- list(TP53 = fit1, RWD = fit2)
    
    g = ggsurvplot(fit,
                   data = tr4,
                   combine = TRUE,
                   palette = "Set1",
                   size = 2,   #line thickness
                   ggtheme = theme_minimal(), #theme
                   font.main = c(15, "black"), #title font
                   font.x = c(15, "bold", "black"), #x-axis font
                   font.y = c(15, "bold", "black"), #y-axis font
                   font.tickslab = c(15, "bold", "black"), #axis numbering font
                   conf.int = FALSE, #confidence interval
                   pval = TRUE, #p-value
                   pval.size = 5, #p-value size
                   pval.coord = c(0, 0.1),
                   risk.table.pos = "out",
                   tables.y.text = FALSE,
                   tables.theme = theme_cleantable(),
                   break.x.by = 2,
                   risk.table.fontsize  = 4,
                   risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                   risk.table.col = "strata", #risk table color
                   risk.table.height = 0.25,
                   ylab = "OS event (%)",
                   xlab = "Time (years)",
                   # xlim = c(0,86),
                   censor = TRUE,
                   censor.shape = 108,
                   censor.size = 4,
                   legend.title = "",
                   legend.labs = c("TP53wt", "TP53mut", "E-RDW (%) low", "E-RDW (%) intermediate", "E-RDW (%) high"),
                   font.legend = c(12, "bold", "black"));g
    
    g = plot_grid(g$plot + guides(colour = guide_legend(nrow = 2)), g$table +  guides(colour = "none"), ncol = 1, align = "v", rel_heights = c(2, 1))
    
    ggsave(plot = g, filename = paste0(results, "/", i, "/KM/RDW_TP53_OS_all_patients.png"), width = 6, height = 6, bg = "white", units = 'in', dpi = 300)
  }
  
  # HDC AML
  if (i == "de novo AML") {
    lab_dg2 = lab_dg1 %>%
      dplyr::inner_join(tr3 %>% dplyr::select(henkilotunnus, HDC, OS_time, OS_event)) %>%
      dplyr::filter(HDC == 1)
    lab_dg2$RDW1 = ntile(lab_dg2$RDW, 3)
    
    # Plot
    fit = survfit(Surv(OS_time, OS_event)~RDW1, data = lab_dg2)
    
    g = ggsurvplot(fit,
                   data = lab_dg2,
                   palette = "Set1",
                   size = 2,   #line thickness
                   ggtheme = theme_minimal(), #theme
                   font.main = c(15, "black"), #title font
                   font.x = c(15, "bold", "black"), #x-axis font
                   font.y = c(15, "bold", "black"), #y-axis font
                   font.tickslab = c(15, "bold", "black"), #axis numbering font
                   conf.int = FALSE, #confidence interval
                   pval = TRUE, #p-value
                   pval.size = 5, #p-value size
                   pval.coord = c(0, 0.1),
                   risk.table.pos = "out",
                   tables.y.text = FALSE,
                   tables.theme = theme_cleantable(),
                   break.x.by = 2,
                   risk.table.fontsize  = 4,
                   risk.table = c("absolute"), #risk table ('absolute', 'percentage', 'abs_pct')
                   risk.table.col = "strata", #risk table color
                   risk.table.height = 0.25,
                   ylab = "OS event (%)",
                   xlab = "Time (years)",
                   # xlim = c(0,86),
                   censor = TRUE,
                   censor.shape = 108,
                   censor.size = 4,
                   font.legend = c(13, "bold", "black"),    #font voi olla esim. "bold" tai "plain"
                   legend.title = "E-RDW (%)",
                   legend.labs = c("Low", "Intermediate", "High")); g
    
    ggsave(plot = g, filename = paste0(results, "/", i, "/KM/RDW_OS_all_patients_HDC.png"), width = 6, height = 6, bg = "white", units = 'in', dpi = 300)
    
  }
  
}


############## TP53 MUTATION TYPE ##############


# Read myelmut
myelmut = readxl::read_xlsx("mounts/research/husdatalake/data/pathology/qpati/myelmut.xlsx") %>%
  dplyr::filter(Gene == "TP53") %>%
  distinct()
## Number of mutations
myelmut_samples = myelmut %>%
  dplyr::arrange(sampletaken) %>%
  dplyr::group_by(henkilotunnus) %>%
  slice(1) %>%
  ungroup()
myelmut1 = myelmut %>%
  dplyr::filter(samplenumber %in% myelmut_samples$samplenumber) %>%
  dplyr::group_by(henkilotunnus, samplenumber) %>%
  mutate(TP53_n = n())
## VAF
myelmut2 = myelmut %>%
  dplyr::filter(samplenumber %in% myelmut_samples$samplenumber) %>%
  dplyr::group_by(henkilotunnus, samplenumber) %>%
  mutate(TP53_vaf = max(VAF, na.rm=TRUE),
         TP53_vaf = gsub(",", ".", TP53_vaf),
         TP53_vaf = as.numeric(TP53_vaf))
## Combine
myelmut3 = full_join(myelmut1 %>% dplyr::select(henkilotunnus, sampletaken, TP53_n),
                     myelmut2 %>% dplyr::select(henkilotunnus, sampletaken, TP53_vaf)) %>%
  distinct()

for (i in c("MDS", "de novo AML")) {
  
  ## Demo
  demo = readRDS(paste0("mounts/research/husdatalake/disease/processed_data/", gsub("de novo ", "", i), "/demographics.rds"))
  
  ## Lab
  lab = readRDS(paste0("mounts/research/husdatalake/disease/data/", gsub("de novo AML", "AML_CML", i), "/multilab.rds"))
  lab1 = lab %>%
    dplyr::filter(henkilotunnus %in% demo$henkilotunnus) %>%
    dplyr::filter(tutkimus_lyhenne == "E -RDW") %>%
    dplyr::select(henkilotunnus, naytteenotto_hetki, tulos)
  
  ## Join with demo
  myelmut4 = myelmut3 %>%
    dplyr::filter(henkilotunnus %in% demo$henkilotunnus) %>%
    dplyr::left_join(lab1) %>%
    dplyr::mutate(timediff = abs(as.Date(naytteenotto_hetki) - as.Date(sampletaken))) %>%
    dplyr::arrange(timediff) %>%
    group_by(henkilotunnus) %>%
    dplyr::slice(1) %>%
    ungroup() %>%
    dplyr::mutate(tulos = as.numeric(tulos),
                  TP53_n = factor(TP53_n))
  
  
  # Plots
  ## Kruskal-Wallis RDW ~ VAF
  g = ggplot(data = myelmut4, aes(x = TP53_n, y = tulos, group=TP53_n)) +
    geom_jitter(width = 0.2, size = 0.5, color = "black") +
    geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
    labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
         y="E-RDW (%)", x="#TP53 mutations") +
    theme_bw() +
    stat_compare_means(method = "kruskal.test", label.x = ifelse(length(levels(myelmut4$TP53_n)) == 3, 1.5,
                                                                 ifelse(length(levels(myelmut4$TP53_n)) == 4, 2, 1.2))) +#, aes(label = ..p.signif..)) +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, colour = "black"),
          plot.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none",
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank()); g
  ggsave(plot = g,
         filename = paste0(results, "/", i, "/RDW_TP53_mutation_number.png"),
         height = 4, width = 4, units = "in", dpi = 300)
  
  ## Correlation RDW ~ VAF
  g = ggplot(data = myelmut4, aes(x = TP53_vaf, y = tulos)) +
    geom_point() +
    labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
         y="E-RDW (%)", x="TP53 VAF (%)") +
    theme_bw() +
    stat_cor(method = "spearman") +
    theme(axis.text.x = element_text(size=12, colour = "black"),
          axis.text.y = element_text(size=12, colour = "black"),
          axis.title=element_text(size=14, colour = "black"),
          plot.title=element_text(size=14, face="bold", colour = "black"),
          legend.position = "none",
          panel.border = element_blank(),
          axis.line = element_line(colour = "black"),
          panel.background = element_blank(),
          panel.grid.minor = element_blank()); g
  ggsave(plot = g,
         filename = paste0(results, "/", i, "/RDW_TP53_mutation_vaf.png"),
         height = 4, width = 4, units = "in", dpi = 300)
}
