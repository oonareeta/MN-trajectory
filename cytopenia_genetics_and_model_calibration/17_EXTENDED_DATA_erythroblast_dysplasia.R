# Calculate missing data

# Libraries
source("mounts/research/src/Rfunctions/library.R")

# General parameters
source("mounts/research/husdatalake/disease/scripts/Preleukemia/parameters")

# Pathology
aml = readRDS("mounts/research/husdatalake/disease/processed_data/AML/demographics.rds") %>%
  dplyr::filter(referral_outside_of_uusimaa_for_allo == FALSE & referral_outside_of_uusimaa == FALSE) %>%
  dplyr::filter(lubridate::year(dg_date_combined) > 2006) %>%
  dplyr::mutate(disease="AML") %>%
  dplyr::select(henkilotunnus, dg_date_combined, disease)
mds = readRDS("mounts/research/husdatalake/disease/processed_data/MDS/demographics.rds") %>%
  dplyr::filter(referral_outside_of_uusimaa_for_allo == FALSE & referral_outside_of_uusimaa == FALSE) %>%
  dplyr::filter(lubridate::year(dg_date_combined) > 2006) %>%
  dplyr::mutate(disease="MDS") %>%
  dplyr::select(henkilotunnus, dg_date_combined, disease)
mf = readRDS("mounts/research/husdatalake/disease/processed_data/MF/demographics.rds") %>%
  dplyr::filter(lubridate::year(dg_date_combined) > 2006) %>%
  dplyr::mutate(disease="MF") %>%
  dplyr::select(henkilotunnus, dg_date_combined, disease)
mn = full_join(full_join(aml, mds), mf)
pat1 = arrow::read_parquet("mounts/research/husdatalake/data/pathology/qpati/pathology_mgg.parquet") %>%
  dplyr::inner_join(mn)
pat = pat1 %>%
  dplyr::mutate(time_diff = abs(as.Date(dg_date_combined) - as.Date(sampletaken))) %>%
  dplyr::arrange(time_diff) %>%
  dplyr::filter(time_diff < 30) %>%
  dplyr::filter(!(statement == "" | is.na(statement))) %>%
  dplyr::filter(!str_detect(statement, "(?i)lisälausunto")) %>%
  dplyr::filter(str_detect(statement, "(?i)(erytro|erythro)")) %>%
  group_by(henkilotunnus) %>%
  slice(1) %>%
  ungroup()

# Extract sentences containing the string "ery"
sentences_with_ery <- str_extract_all(pat$statement, "\\b[^.!?]*?ery[^.!?]*[.!?]")

# Unlist and get unique sentences
unique_sentences <- unique(unlist(sentences_with_ery))

# Store the unique sentences in a list
sentence_list <- as.list(unique_sentences)


# Function to extract sentences containing specified keywords
extract_sentences <- function(text, keywords) {
  # Split the text into sentences
  sentences <- str_split(text, "\\.")[[1]]
  
  # Filter sentences containing keywords
  selected_sentences <- sentences[str_detect(sentences, paste(keywords, collapse = "|", sep = ""))]
  
  # Return the selected sentences
  return(selected_sentences)
}

# Specify keywords
keywords0 = "(?i)ery"
keywords1 <- c("(?i)(silmikointi|deform|monituma|dysplas|kromatiini|megalo|kypsymishä|tumamuut|vakuol|patolog)")
keywords2 = c("(?i)(((ring|sidero)[^.]*(yli |[[:digit:]]{2}%))|((yli |([[:digit:]]{2}%))[^.]*(ring|sidero)))")
keywords3 <- c("(?i)(ei [^.]{0,20}(viittaavaa|todeta)|pienessä|yksitt|lievä|lieviä|vähäi|ani[^.]{0,1}harva| ei [^.]{0,20}sel| ei [^.]{0,20}vaku| ei [^.]{0,20}tode|(ring|sidero)[^.]*(^(alle|<))|(^(alle|<))[^.]*(ring|sidero))")

# Apply the function to each row of the dataframe
result1 <- pat %>%
  rowwise() %>%
  mutate(selected_sentences1 = list(extract_sentences(statement, keywords0))) %>%
  unnest(selected_sentences1)

result1 <- result1 %>%
  dplyr::filter(str_detect(selected_sentences1, keywords1) | str_detect(selected_sentences1, keywords2))

result2 <- result1 %>%
  dplyr::filter(!(str_detect(selected_sentences1, keywords3) & !str_detect(selected_sentences1, keywords2)))

# Final
ery_final = pat %>% 
  dplyr::mutate(ery_dys = ifelse(henkilotunnus %in% result2$henkilotunnus, TRUE, FALSE)) %>%
  dplyr::distinct(henkilotunnus, ery_dys)

# Save
saveRDS(ery_final, paste0(export, "/Erythroid_dysplasia.rds"))


# Dysplasia
mgg1 = readxl::read_xlsx("mounts/research/husdatalake/data/pathology/qpati/dysplasia.xlsx") %>%
  dplyr::filter(henkilotunnus %in% mn$henkilotunnus) %>%
  dplyr::select(henkilotunnus, samplenumber, sampletaken,
                starts_with("ERY"), starts_with("APL"), starts_with("MEG"), 
                MGK_count, cellularity, bm_blast, bm_ly, bm_plasma, bm_sideroblast) %>%
  dplyr::left_join(mn) %>%
  mutate(time_diff = abs(as.Date(sampletaken) - as.Date(dg_date_combined))) %>%
  arrange(time_diff) %>%
  group_by(henkilotunnus) %>%
  dplyr::slice(1) %>%
  ungroup()
## Fill NA with 0
mgg1 <- mgg1 %>%
  mutate(across(starts_with(c("ERY", "APL", "MEG")), ~ replace_na(., 0)))
## Process
mgg1 = mgg1 %>%
  dplyr::mutate(ERY_vacuolated_mild = ifelse(ERY_vacuolated_mild == 1, 0.5, ERY_vacuolated),
                ERY_large_ery_mild = ifelse(ERY_large_ery_mild == 1, 0.5, ERY_large_ery),
                ERY_megaloblast_mild = ifelse(ERY_megaloblast_mild == 1, 0.5, ERY_megaloblast),
                ERY_dysmorphic_mild = ifelse(ERY_dysmorphic_mild == 1, 0.5, ERY_dysmorphic),
                ERY_multinucleated_mild = ifelse(ERY_multinucleated_mild == 1, 0.5, ERY_multinucleated),
                ERY_dysplasia = pmax(ERY_vacuolated_mild, ERY_large_ery_mild, ERY_megaloblast_mild, ERY_dysmorphic_mild, ERY_multinucleated_mild, na.rm=TRUE),
                APL_dysplasia = pmax(APL_auer, APL_nucleus, APL_hypergran, na.rm=TRUE),
                MGK_dysplasia = pmax(MEG_separated, MEG_micro, MEG_large, MEG_hyperlob, MEG_hypolob, na.rm=TRUE),
                bm_blast = ifelse(is.na(bm_blast) | bm_blast == "Norm", 2.5, bm_blast),
                bm_ly = ifelse(is.na(bm_ly) | bm_ly == "Norm", 30, bm_ly),
                bm_plasma = ifelse(is.na(bm_plasma) | bm_plasma == "Norm", 2.5, bm_plasma),
                bm_sideroblast = ifelse(is.na(bm_sideroblast) | bm_sideroblast == "Norm", 0, bm_sideroblast))




# Plot
for (i in c("MDS", "de novo AML", "MF")) {
  
  print(i)
  
  ## Demo
  demo = readRDS(paste0("mounts/research/husdatalake/disease/processed_data/", gsub("de novo ", "", i), "/demographics.rds"))
  
  ## Lab
  lab = readRDS(paste0("mounts/research/husdatalake/disease/processed_data/", gsub("de novo ", "", i), "/lab_dg.rds"))
  lab1 = lab %>%
    dplyr::filter(henkilotunnus %in% demo$henkilotunnus) %>%
    dplyr::filter(tutkimus_lyhenne == "E -RDW") %>%
    dplyr::select(henkilotunnus, tulos=tulos_mean)
  lab1 = lab1 %>%
    dplyr::left_join(ery_final) %>%
    dplyr::filter(!is.na(ery_dys))
  
  g = ggplot(data = lab1, aes(x = ery_dys, y = tulos)) + #, group=time_to_dg_mo
    geom_jitter(width = 0.2, size = 0.5, color = "black") +
    geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
    labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
         y="E-RDW (%)", x="Erythroid dysplasia") +
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
         filename = paste0(results, "/", i, "/RDW_erythroid_dysplasia_", i, ".png"),
         height = 4, width = 4, units = "in", dpi = 300)
}


# Plot
for (i in c("MDS", "de novo AML", "MF")) {
  
  print(i)
  
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  # Read lab data
  df = fread(paste0(export, "/lab_data_for_modelling_", i, ".csv")) %>%
    dplyr::filter(henkilotunnus %in% unique(df1$henkilotunnus)) %>%
    dplyr::filter(time_to_dg<-90 & time_to_dg>-730)
  
  # Keep RDW
  lab1 = df %>%
    dplyr::filter((tutkimus_lyhenne_yksikko == "E-RDW (%)")) %>%
    # Summarise median by patient and laboratory test
    dplyr::group_by(henkilotunnus) %>%
    mutate(tulos = mean(tulos, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(henkilotunnus, tulos)
  
  lab1 = lab1 %>%
    dplyr::left_join(ery_final) %>%
    dplyr::filter(!is.na(ery_dys))
  
  g = ggplot(data = lab1, aes(x = ery_dys, y = tulos)) + #, group=time_to_dg_mo
    geom_jitter(width = 0.2, size = 0.5, color = "black") +
    geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
    labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
         y="E-RDW (%)", x="Erythroid dysplasia") +
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
         filename = paste0(results, "/", i, "/RDW_erythroid_dysplasia_", i, "_preleukemia_labs.png"),
         height = 4, width = 4, units = "in", dpi = 300)
}


# The same with the exact dysplasia
# Plot
for (i in c("MDS", "de novo AML", "MF")) {
  
  print(i)
  
  ## Demo
  demo = readRDS(paste0("mounts/research/husdatalake/disease/processed_data/", gsub("de novo ", "", i), "/demographics.rds"))
  
  ## Lab
  lab = readRDS(paste0("mounts/research/husdatalake/disease/processed_data/", gsub("de novo ", "", i), "/lab_dg.rds"))
  lab1 = lab %>%
    dplyr::filter(henkilotunnus %in% demo$henkilotunnus) %>%
    dplyr::filter(tutkimus_lyhenne == "E -RDW") %>%
    dplyr::select(henkilotunnus, tulos=tulos_mean)
  lab1 = lab1 %>%
    dplyr::left_join(mgg1)
  
  for (j in c("ERY_vacuolated_mild", "ERY_large_ery_mild", "ERY_megaloblast_mild", "ERY_dysmorphic_mild", "ERY_multinucleated_mild")) {
    
    lab1$tmp = lab1[[j]]
    j1 = to_sentence_case(gsub("ERY_|_mild|_ery", "", j))
    
    if (lab1 %>%
        dplyr::filter(!is.na(tmp)) %>%
        distinct(tmp) %>%
        nrow() > 1 & lab1 %>%
        dplyr::filter(!is.na(tmp)) %>%
        dplyr::filter(tmp > 0.5) %>%
        distinct(tmp) %>%
        nrow() > 1) {
      
      g = ggplot(data = lab1 %>%
                   dplyr::filter(!is.na(tmp)) %>%
                   mutate(tmp = ifelse(tmp > 0.5, TRUE, FALSE),
                          tmp = factor(tmp)), aes(x = tmp, y = tulos)) + #, group=time_to_dg_mo
        geom_jitter(width = 0.2, size = 0.5, color = "black") +
        geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
        labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
             y="E-RDW (%)", x=j1) +
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
             filename = paste0(results, "/", i, "/RDW_erythroid_dysplasia_", i, "_", j1, ".png"),
             height = 4, width = 4, units = "in", dpi = 300)
    }
    
  }
}


# Plot
for (i in c("MDS", "de novo AML", "MF")) {
  
  print(i)
  
  df1 = readRDS(paste0(export, "/", i, "_lab_demo.rds"))
  
  # Read lab data
  df = fread(paste0(export, "/lab_data_for_modelling_", i, ".csv")) %>%
    dplyr::filter(henkilotunnus %in% unique(df1$henkilotunnus)) %>%
    dplyr::filter(time_to_dg<-90 & time_to_dg>-730)
  
  # Keep RDW
  lab1 = df %>%
    dplyr::filter((tutkimus_lyhenne_yksikko == "E-RDW (%)")) %>%
    # Summarise median by patient and laboratory test
    dplyr::group_by(henkilotunnus) %>%
    mutate(tulos = mean(tulos, na.rm=TRUE)) %>%
    slice(1) %>%
    ungroup() %>%
    dplyr::select(henkilotunnus, tulos)
  
  lab1 = lab1 %>%
    dplyr::left_join(mgg1)
  
  for (j in c("ERY_vacuolated_mild", "ERY_large_ery_mild", "ERY_megaloblast_mild", "ERY_dysmorphic_mild", "ERY_multinucleated_mild")) {
    
    lab1$tmp = lab1[[j]]
    j1 = to_sentence_case(gsub("ERY_|_mild|_ery", "", j))
    
    if (lab1 %>%
        dplyr::filter(!is.na(tmp)) %>%
        distinct(tmp) %>%
        nrow() > 1) {
      
      g = ggplot(data = lab1 %>%
                   dplyr::filter(!is.na(tmp)) %>%
                   mutate(tmp = ifelse(tmp > 0.5, TRUE, FALSE),
                          tmp = factor(tmp)), aes(x = tmp, y = tulos)) + #, group=time_to_dg_mo
        geom_jitter(width = 0.2, size = 0.5, color = "black") +
        geom_boxplot(alpha = 0.3, color = "#e41a1c", outliers = FALSE) +
        labs(title = gsub("de novo AML", "De novo AML", gsub("MF", "Primary MF", i)),
             y="E-RDW (%)", x=j1) +
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
             filename = paste0(results, "/", i, "/RDW_erythroid_dysplasia_", i, "_" , j1, "_preleukemia_labs.png"),
             height = 4, width = 4, units = "in", dpi = 300)
    }
  }
  
}
