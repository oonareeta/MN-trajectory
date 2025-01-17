# Calculate missing data


# Author: Oscar Brück

# Libraries
source("./library.R")


# General parameters
source("./parameters")


############################################################


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
pat1 = readRDS("mounts/research/husdatalake/data/pathology/qpati/pathology_mgg.rds") %>%
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
  
  # Exclude patients with wrong diagnosis
  exclude = readxl::read_xlsx(paste0("mounts/research/husdatalake/disease/general/exclude_pts.xlsx"))
  if (i == "de novo AML") {
    j = "AML"
  } else {
    j = i
  }
  df1 = df1 %>%
    dplyr::filter(!(disease == "AML" & henkilotunnus %in% exclude[exclude$disease=="AML",]$henkilotunnus))
  df1 = df1 %>%
    dplyr::filter(!(disease == "MF" & henkilotunnus %in% exclude[exclude$disease=="MF",]$henkilotunnus))
  df1 = df1 %>%
    dplyr::filter(!(disease == "MDS" & henkilotunnus %in% exclude[exclude$disease=="MDS",]$henkilotunnus))
  
  # Remove secondary MF
  pmf = readRDS("mounts/research/husdatalake/disease/processed_data/MF/secondary_MF.rds") %>%
    dplyr::filter(!disease == "Primary")
  df1 = df1 %>%
    dplyr::filter(!(disease == "MF" & henkilotunnus %in% pmf$henkilotunnus))
  
  
  # Read lab data
  df = fread(paste0(export, "/lab_data_for_modelling_", i, ".csv")) %>%
    dplyr::filter(henkilotunnus %in% unique(df1$henkilotunnus)) %>%
    dplyr::filter(time_to_dg<-30 & time_to_dg>-730)
  
  
  # Remove PV and ET patients identified based on JAK2 mutations and diagnosis table
  df = df %>%
    dplyr::filter(!henkilotunnus %in% c('02139_43396117','02139_84320876','02139_05958098','02139_09649728','02139_43105799','02139_86518942','02139_58286011',
                                        '02139_74536441','02139_04048905','02139_66454572','02139_84083902','02139_08144748','02139_93742663','02139_09325470',
                                        '02139_98606703','02139_32936070','02139_74525166','02139_01261095','02139_72491838','02139_16511870','02139_04981848',
                                        '02139_72056477','02139_47431683','02139_91416073'))
  
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
