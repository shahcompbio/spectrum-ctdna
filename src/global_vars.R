## plotting themes --------------------------------

theme_cowplot2 <- function(...) {
  theme_cowplot(font_size = 16, font_family = "sans", ...) %+replace%
    theme(strip.background = element_blank(),
          panel.background = element_rect(fill = "transparent", color = NA),
          plot.background = element_rect(fill = "transparent", color = NA),
          panel.border = element_blank())
}
theme_set(theme_cowplot2())

remove_xaxis <- theme(axis.title.x = element_blank(),
                      axis.text.x = element_blank(),
                      axis.ticks.x = element_blank(),
                      axis.line.x = element_blank())

remove_yaxis <- theme(axis.title.y = element_blank(),
                      axis.text.y = element_blank(),
                      axis.ticks.y = element_blank(),
                      axis.line.y = element_blank())

remove_guides <- guides(color = F, fill = F, shape = F, alpha = F, size = F)


## ggsave wrapper suppressing dingbats symbols 
## for adobe illustrator compatibility
ggsave_pdf <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = "in",#units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, ...) {
  ggsave(filename = filename, plot = plot, device = cairo_pdf, path = path, 
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, ...)
}

ggsave_png <- function(filename, plot = last_plot(), device = NULL, path = NULL, 
                       scale = 1, width = NA, height = NA, units = "in",#units = c("in", "cm", "mm"), 
                       dpi = 300, limitsize = TRUE, type = "cairo", ...) {
  ggsave(filename = filename, plot = plot, device = device, path = path, 
         scale = scale, width = width, height = height, units = units, 
         dpi = dpi, limitsize = limitsize, type = type, ...)
}

## load color code --------------------------------

clrs <- yaml::read_yaml("resources/annotation/colors.yaml") %>%
  lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

clrs$patient_id_short <- clrs$patient_id
names(clrs$patient_id_short) <- str_remove_all(names(clrs$patient_id), "SPECTRUM-OV-")

shps <- yaml::read_yaml("resources/annotation/shapes.yaml") %>% 
  lapply(function(x) map_depth(x, vec_depth(x)-2, unlist))

## load database ----------------------------------

db <- readr::read_rds("resources/db/ctdna/SPECTRUM.rds")

## define patients included in the study -----------

# Define confirmed HGS patients
hgsoc_patients <- db$gyn_diagnosis %>%
  filter(gyn_diagnosis_histology == "HGS") %>%
  pull(patient_id)

# Define patients on clinical trials
protocol_patients <- db$consents %>%
  filter(patient_consent_irb == "17-182") %>%
  pull(patient_id)

# Create list of patients on the SPECTRUM TME study
# - Exclude non-HGSOC patients
# - Exclude patients on clinical trials (e.g. 17-182)
included_patients <- db$patients %>%
  filter(patient_inclusion_exclusion=="Included") %>%
  # filter(patient_cohort_version___2=="Checked") %>%
  # filter(patient_id %in% hgsoc_patients) %>%
  # filter(!patient_id %in% protocol_patients) %>%
  pull(patient_id)

# Define patients included in the study with scDNA data
cfdna_patients <- db$sequencing_cfdna %>%
  filter(patient_id %in% included_patients) %>%
         # qc_status == "Pass") %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with scDNA data
scdna_patients <- db$sequencing_scdna %>%
  filter(patient_id %in% included_patients,
         qc_status == "Pass") %>%
  pull(patient_id) %>%
  unique

# # Define patients included in the study with scRNA data
# scrna_patients <- db$sequencing_scrna %>%
#   # filter(patient_id %in% included_patients) %>%
#   filter(qc_status == "Pass") %>%
#   pull(patient_id) %>%
#   unique

# Define patients included in the study with bulk DNA data
bulk_dna_patients <- db$sequencing_bulk_dna %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scdna_patients, cfdna_patients),
         qc_status == "Pass") %>%
  distinct(patient_id) %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with Nanopore data
nanopore_patients <- db$sequencing_nanopore %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scdna_patients, cfdna_patients),
         qc_status == "Pass") %>%
  distinct(patient_id) %>%
  pull(patient_id) %>%
  unique

# Define patients included in the study with IMPACT data
impact_patients <- db$sequencing_msk_impact_custom %>%
  filter(patient_id %in% included_patients,
         patient_id %in% union(scdna_patients, cfdna_patients),
         qc_status == "Pass") %>%
  distinct(patient_id) %>%
  pull(patient_id) %>%
  unique

## load mutational signatures ----------------------

signature_tbl <- db$mutational_signatures %>%
  mutate(consensus_signature = ordered(consensus_signature, levels = names(clrs$consensus_signature))) %>% 
  arrange(patient_id)

## load scDNA meta data -----------------------------

scdna_meta_tbl <- db$sequencing_scdna %>% 
  filter(patient_id %in% included_patients) %>% 
  dplyr::rename(sample = isabl_id) %>% 
  distinct(sample, .keep_all = T) %>% 
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>%
  left_join(signature_tbl, by = "patient_id")

## load cfDNA meta data -------------------------------

cfdna_meta_tbl <- db$sequencing_cfdna %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-")) %>% 
  filter(patient_id %in% included_patients) %>% 
  # left_join(db$specimens_lab_medicine %>%
  #             mutate(specimen_lab_medicine_accession = as.character(specimen_lab_medicine_accession)) %>%
  #             select(-c("patient_id")), 
  #           by = c("isabl_id" = "specimen_lab_medicine_accession")) %>%
  left_join(signature_tbl, by = "patient_id")

## load WGS meta data -------------------------------

bulk_dna_meta_tbl <- db$sequencing_bulk_dna %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         sample_id_short = str_remove_all(sample_id, "OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  filter(patient_id %in% included_patients) %>% 
  left_join(signature_tbl, by = "patient_id")

## load MSK-IMPACT meta data ------------------------

impact_meta_tbl <- db$sequencing_msk_impact_custom %>%
  mutate(patient_id_short = str_remove_all(patient_id, "SPECTRUM-OV-"),
         tumor_supersite = str_replace_all(tumor_supersite, "Upper Quadrant", "UQ")) %>% 
  mutate(tumor_megasite = ifelse(!tumor_supersite %in% c("Adnexa", "Ascites"),
                                 "Other", tumor_supersite)) %>% 
  mutate(tumor_supersite = ordered(tumor_supersite, levels = names(clrs$tumor_supersite))) %>% 
  filter(patient_id %in% included_patients) %>% 
  left_join(signature_tbl, by = "patient_id")
