# Import external required packages
library(BiocManager)
library(cancereffectsizeR)
library(data.table)
library(tidyverse)
library(dplyr)
library(readr)
# library(maftools)
library(ggpmisc)
library(openxlsx)
library(readxl)
library(survminer)
library(survival)
library("maxstat")
library(BSgenome.Hsapiens.UCSC.hg19)
# library(SigProfilerAssignmentR)
library(cowplot)
# library(deconstructSigs)
library(cutpointr)
library(gtools)
library(scales)
library(grid)
library(gridExtra)
library(contsurvplot)
library(pammtools)
library(riskRegression)

# setwd("files")

maf_format <- c("Hugo_Symbol", "Chromosome", "Start_Position", "End_Position", "Reference_Allele", "Variant_Classification", "Tumor_Sample_Barcode")


#----------Import data from George et al., 2015 (Nature)------------------


George_et_al <- "george_et_al.xlsx"

if (!file.exists(George_et_al)) {
  george_et_al_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fnature14664/MediaObjects/41586_2015_BFnature14664_MOESM72_ESM.xlsx"
  download.file(george_et_al_url, destfile = George_et_al, mode = "wb")
}


## WXS data
george_et_al <- read.xlsx(George_et_al, sheet = 3, startRow = 3, cols = 1:23)
length(unique(george_et_al$PAT_ID))

setnames(george_et_al, old = "Start", new = "Start_Position", skip_absent=TRUE)
setnames(george_et_al, old = "End", new = "End_Position", skip_absent=TRUE)
setnames(george_et_al, old = "PAT_ID", new = "Unique_Patient_Identifier", skip_absent=TRUE)
setnames(george_et_al, old = "Wild_Type", new = "Reference_Allele", skip_absent=TRUE)
setnames(george_et_al, old = "Mutant", new = "Tumor_Seq_Allele", skip_absent=TRUE)
setnames(george_et_al, old = "Gene_Hugo", new = "Hugo_Symbol", skip_absent=TRUE)
setnames(george_et_al, old = "Type_1", new = "Variant_Classification", skip_absent=TRUE)
setnames(george_et_al, old = "Mut_ID", new = "Tumor_Sample_Barcode", skip_absent=TRUE)


george_et_al_maf <- preload_maf(maf = george_et_al, refset = "ces.refset.hg19", tumor_allele_col = "Tumor_Seq_Allele", keep_extra_columns = maf_format)

## Clinical data
george_et_al_clinical <- read.xlsx(George_et_al, sheet = 1, startRow = 4, cols = c(1, 12:29))

# Remove unwanted columns
keep_cols <- c("Sample-ID", "age", "sex", "stage_UICC", "smoking_status", "previous.therapeutic.treatment.for.SCLC", "chemotherapy.(yes/no)","Status.(at.time.of.last.follow-up)", "progression-free_survival.(months)", "overall_survival.(months)")
george_et_al_clinical <- george_et_al_clinical[, keep_cols]

# Rename columns
setnames(george_et_al_clinical, old = "Sample-ID", new = "Unique_Patient_Identifier", skip_absent=TRUE)
setnames(george_et_al_clinical, old = "stage_UICC", new = "staging", skip_absent=TRUE)


george_et_al_clinical <- george_et_al_clinical[
  george_et_al_clinical$Unique_Patient_Identifier %in% unique(george_et_al_maf$Unique_Patient_Identifier), ]

george_et_al_clinical <- george_et_al_clinical %>%
  mutate(
    previous.therapeutic.treatment.for.SCLC = recode(previous.therapeutic.treatment.for.SCLC, "untreated"= FALSE, "relapse" = TRUE),
    smoking_status = if_else(is.na(smoking_status), "non-smoker", "smoker")
  )


#--------Import data from Rudin et al., 2012 (Nature)---------------------


Rudin_et_al <- "rudin_et_al.xls"

if (!file.exists(Rudin_et_al)) {
  rudin_et_al_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fng.2405/MediaObjects/41588_2012_BFng2405_MOESM20_ESM.xls"
  download.file(rudin_et_al_url, destfile = Rudin_et_al, mode = "wb")
}

## WXS Data
rudin_et_al <- read_excel(Rudin_et_al, sheet = 4, skip = 1)
rudin_et_al <- rudin_et_al[, 1:18]  # Keep only relevant columns
rudin_et_al$End_Position <- NA


setnames(rudin_et_al, old = "Position", new = "Start_Position", skip_absent=TRUE)
setnames(rudin_et_al, old = "Tumor_Normal Sample ID", new = "Unique_Patient_Identifier", skip_absent=TRUE)
rudin_et_al$Tumor_Sample_Barcode <- rudin_et_al$Unique_Patient_Identifier
setnames(rudin_et_al, old = "Ref", new = "Reference_Allele", skip_absent=TRUE)
setnames(rudin_et_al, old = "Var", new = "Tumor_Seq_Allele", skip_absent=TRUE)
setnames(rudin_et_al, old = "GeneName", new = "Hugo_Symbol", skip_absent=TRUE)
setnames(rudin_et_al, old = "Sift", new = "Variant_Classification", skip_absent=TRUE)



rudin_et_al_maf <- preload_maf(maf = rudin_et_al, refset = "ces.refset.hg19", tumor_allele_col = "Tumor_Seq_Allele", keep_extra_columns = maf_format)

## Clinical Data Unavailable


# # -----------------Import data from Msk_Met_2021--------------------------
# 
# 
# if (!file.exists("msk_met_2021.tar.gz")) {
#   msk_met_url <- "https://cbioportal-datahub.s3.amazonaws.com/msk_met_2021.tar.gz"
#   download.file(msk_met_url, destfile = "msk_met_2021.tar.gz", mode = "wb")
# }
# 
# untar("msk_met_2021.tar.gz", exdir = ".")
# 
# ## Targeted Sequencing 
# msk_met_all_txt_file <- read_tsv("msk_met_2021/data_mutations.txt", show_col_types = FALSE) # Not specific to SCLC
# 
# # colnames(msk_met_all_txt_file)
# msk_met_all_txt_file <- msk_met_all_txt_file %>%
#   mutate(Unique_Patient_Identifier = str_extract(Tumor_Sample_Barcode, "^[^-]+-[^-]+"))
# 
# # msk_met_all_maf <- preload_maf(maf = msk_met_all_txt_file, refset = "ces.refset.hg19") # Not specific to SCLC
# 
# ## Clinical Data
# 
# # View files for downloaded msk_met clinical data
# 
# # Set file path
# file_path <- "Msk_met_all_clinical.xlsx" # NIHMS1774757-supplement-5.xlsx from Nguyen et al., 2022 (Cell)
# 
# # Read the Excel file
# msk_met_all_clinical <- read.xlsx(file_path, sheet = 2, startRow = 3, cols = 1:28)
# 
# # Filter for Small Cell Lung Cancer
# msk_met_clinical_txt_file <- msk_met_all_clinical[msk_met_all_clinical$cancer_type == "small cell lung cancer", ]
# 
# # Remove unwanted columns
# keep_cols <- c("patient_id", "prior_treatement", "seq_report_age", "os_days", "os_status", "sex")
# msk_met_clinical_txt_file <- msk_met_clinical_txt_file[, keep_cols]
# 
# # Rename columns
# names(msk_met_clinical_txt_file) <- c("Unique_Patient_Identifier", 
#                                       "previous.therapeutic.treatment.for.SCLC", 
#                                       "age", 
#                                       "overall_survival.(months)", 
#                                       "Status.(at.time.of.last.follow-up)",
#                                       "sex")
# 
# # Extract patient IDs
# sclc_patient_ids <- msk_met_clinical_txt_file[["Unique_Patient_Identifier"]]
# 
# # Convert units / clean columns
# msk_met_clinical_txt_file$`overall_survival.(months)` <-  msk_met_clinical_txt_file$`overall_survival.(months)` / 30
# msk_met_clinical_txt_file$age <- msk_met_clinical_txt_file$age / 365
# msk_met_clinical_txt_file$sex <- tolower(msk_met_clinical_txt_file$sex)
# 
# ## Filter msk_met_maf using clinical data to be SCLC specific
# msk_met_clinical_txt_file <- msk_met_clinical_txt_file[msk_met_clinical_txt_file$Unique_Patient_Identifier != "P-0032529", ] # removing for sigprofilerassignmentr analysis
# sclc_patient_ids <- msk_met_clinical_txt_file$Unique_Patient_Identifier
# nrow(msk_met_clinical_txt_file)
# msk_met_txt_file <- msk_met_all_txt_file %>%
#   filter(Unique_Patient_Identifier %in% sclc_patient_ids)
# # We will use the following maf file to guide clinical data
# # msk_met_maf <- preload_maf(maf = msk_met_txt_file, refset = "ces.refset.hg19")
# 
# # Define coverage for targeted data set
# 
#   # Import the classes of coverage used by msk_met_2021 
# IMPACT468 <- readLines("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/reference_data/gene_panels/data_gene_panel_impact468.txt")
# # Isolate the row that starts with "gene_list:" as a character
# IMPACT468_genes <- IMPACT468[grepl("^gene_list:", IMPACT468)]
# # Delete the "gene_list:" from the character
# IMPACT468_genes <- sub("^gene_list:\\s*", "", IMPACT468_genes)
# # Replace the variable gaps between gene names with a consistent singular space
# IMPACT468_genes <- gsub("\\s+", " ", IMPACT468_genes)
# # Use the consistent singular space gap to separate gene names
# IMPACT468_genes <- unlist(strsplit(trimws(IMPACT468_genes), " "))
# 
# IMPACT341 <- readLines("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/reference_data/gene_panels/data_gene_panel_impact341.txt")
# # Isolate the row that starts with "gene_list:" as a character
# IMPACT341_genes <- IMPACT341[grepl("^gene_list:", IMPACT341)]
# # Delete the "gene_list:" from the character
# IMPACT341_genes <- sub("^gene_list:\\s*", "", IMPACT341_genes)
# # Replace the variable gaps between gene names with a consistent singular space
# IMPACT341_genes <- gsub("\\s+", " ", IMPACT341_genes)
# # Use the consistent singular space gap to separate gene names
# IMPACT341_genes <- unlist(strsplit(trimws(IMPACT341_genes), " "))
# 
# IMPACT410 <- readLines("https://media.githubusercontent.com/media/cBioPortal/datahub/refs/heads/master/reference_data/gene_panels/data_gene_panel_impact410.txt")
# # Isolate the row that starts with "gene_list:" as a character
# IMPACT410_genes <- IMPACT410[grepl("^gene_list:", IMPACT410)]
# # Delete the "gene_list:" from the character
# IMPACT410_genes <- sub("^gene_list:\\s*", "", IMPACT410_genes)
# # Replace the variable gaps between gene names with a consistent singular space
# IMPACT410_genes <- gsub("\\s+", " ", IMPACT410_genes)
# # Use the consistent singular space gap to separate gene names
# IMPACT410_genes <- unlist(strsplit(trimws(IMPACT410_genes), " "))
# 
# # Import a data.frame that includes the specification of the coverage class used for each patient
# msk_met_gene_panel_all_matrix <- read_tsv("msk_met_2021/data_gene_panel_matrix.txt", show_col_types = FALSE) # Not specific to SCLC
# 
# keep_cols <- c("SAMPLE_ID", "mutations")
# msk_met_gene_panel_all_matrix <- msk_met_gene_panel_all_matrix[, keep_cols]
# 
# names(msk_met_gene_panel_all_matrix) <- c("Tumor_Sample_Barcode", 
#                                 "source")
# msk_met_gene_panel_all_matrix$Unique_Patient_Identifier <- sub("^(([^-]+-[^-]+)).*", "\\1", msk_met_gene_panel_all_matrix$Tumor_Sample_Barcode)
# # Filter to only include SCLC patients
# msk_met_gene_panel_matrix <- msk_met_gene_panel_all_matrix[msk_met_gene_panel_all_matrix$Unique_Patient_Identifier %in% sclc_patient_ids, ]
# 
# # Use the specification of coverage class to filter msk_met data and then preload into maf format
# samples_IMPACT468 <- msk_met_gene_panel_matrix %>%
#   filter(source == "IMPACT468") %>%
#   pull(Tumor_Sample_Barcode)
# IMPACT468_msk_met_txt_file <- msk_met_txt_file %>%
#   filter(Tumor_Sample_Barcode %in% samples_IMPACT468)
# IMPACT468_msk_met_maf <- preload_maf(maf = IMPACT468_msk_met_txt_file, refset = "ces.refset.hg19", keep_extra_columns = maf_format)
# IMPACT468_msk_met_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$names %in% IMPACT468_genes]
# 
# samples_IMPACT341 <- msk_met_gene_panel_matrix %>%
#   filter(source == "IMPACT341") %>%
#   pull(Tumor_Sample_Barcode)
# IMPACT341_msk_met_txt_file <- msk_met_txt_file %>%
#   filter(Tumor_Sample_Barcode %in% samples_IMPACT341)
# IMPACT341_msk_met_maf <- preload_maf(maf = IMPACT341_msk_met_txt_file, refset = "ces.refset.hg19", keep_extra_columns = maf_format)
# IMPACT341_msk_met_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$names %in% IMPACT341_genes]
# 
# samples_IMPACT410 <- msk_met_gene_panel_matrix %>%
#   filter(source == "IMPACT410") %>%
#   pull(Tumor_Sample_Barcode)
# IMPACT410_msk_met_txt_file <- msk_met_txt_file %>%
#   filter(Tumor_Sample_Barcode %in% samples_IMPACT410)
# IMPACT410_msk_met_maf <- preload_maf(maf = IMPACT410_msk_met_txt_file, refset = "ces.refset.hg19", keep_extra_columns = maf_format)
# IMPACT410_msk_met_coverage <- ces.refset.hg19$gr_genes[ces.refset.hg19$gr_genes$names %in% IMPACT410_genes]
# 
# -------------------Import External SCLC Data ---------------------------

## Jiang et al., 2016 (PLOS Genetics)

### WXS
Jiang_et_al_WXS <- "jiang_et_al_WXS.xlsx"

if (!file.exists(Jiang_et_al_WXS)) {
  jiang_et_al_WXS_url <- "https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.1005895.s019"
  download.file(jiang_et_al_WXS_url, destfile = Jiang_et_al_WXS, mode = "wb")
}

jiang_et_al <- read.xlsx(Jiang_et_al_WXS, sheet = 1, startRow = 2, cols = 1:22)

setnames(jiang_et_al, old = "Start", new = "Start_Position", skip_absent=TRUE)
setnames(jiang_et_al, old = "Patient", new = "Unique_Patient_Identifier", skip_absent=TRUE)
setnames(jiang_et_al, old = "Ref", new = "Reference_Allele", skip_absent=TRUE)
setnames(jiang_et_al, old = "Alt", new = "Tumor_Seq_Allele", skip_absent=TRUE)
setnames(jiang_et_al, old = "Chr", new = "Chromosome", skip_absent=TRUE)
setnames(jiang_et_al, old = "Gene.refGene", new = "Hugo_Symbol", skip_absent=TRUE)

jiang_et_al_maf <- preload_maf(maf = jiang_et_al, refset = "ces.refset.hg19", tumor_allele_col = "Tumor_Seq_Allele", keep_extra_columns = "Hugo_Symbol")


### Clinical data
Jiang_et_al_clinical <- "jiang_et_al_clinical.xlsx"

if (!file.exists(Jiang_et_al_clinical)) {
  jiang_et_al_clinical_url <- "https://journals.plos.org/plosgenetics/article/file?type=supplementary&id=10.1371/journal.pgen.1005895.s012"
  download.file(jiang_et_al_clinical_url, destfile = Jiang_et_al_clinical, mode = "wb")
}

jiang_et_al_clinical <- read.xlsx(Jiang_et_al_clinical, sheet = 1, startRow = 2, cols = c(1, 6:17))
colnames(jiang_et_al_clinical)[1] <- "Unique_Patient_Identifier"

# Remove unwanted columns
keep_cols <- c("Unique_Patient_Identifier", "Age", "Gender", "UICC.stage", "Smoke", "Months", "Survival.Status", "chemotherapy-exposed.biopsy")
jiang_et_al_clinical <- jiang_et_al_clinical[, keep_cols]

# Rename Columns
names(jiang_et_al_clinical) <- c("Unique_Patient_Identifier", 
                                 "age", 
                                 "sex", "staging",
                                 "smoking_status",
                                 "overall_survival.(months)", 
                                 "Status.(at.time.of.last.follow-up)", "chemotherapy.(yes/no)")

# Clean up data into appropriate format

jiang_et_al_clinical <- jiang_et_al_clinical %>%
  mutate(
    sex = recode(sex, `1` = "male", `2` = "female"),
    smoking_status = recode(smoking_status, `1` = "smoker", `2` = "non-smoker"),
    `Status.(at.time.of.last.follow-up)` = recode(`Status.(at.time.of.last.follow-up)`, `0` = "dead", `1` = "alive"),
    `chemotherapy.(yes/no)` = recode(`chemotherapy.(yes/no)`, `Chemo` = "yes", `Naive` = "no")
  )

jiang_et_al_clinical <- jiang_et_al_clinical %>%
  filter(Unique_Patient_Identifier %in% unique(jiang_et_al_maf$Unique_Patient_Identifier))


## Zhou et al., 2021 (Nature)

### WXS
Zhou_et_al_WXS <- "zhou_et_al_WXS.xlsx"

if (!file.exists(Zhou_et_al_WXS)) {
  zhou_et_al_WXS_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-25787-x/MediaObjects/41467_2021_25787_MOESM4_ESM.xlsx"
  download.file(zhou_et_al_WXS_url, destfile = Zhou_et_al_WXS, mode = "wb")
}

zhou_et_al <- read.xlsx(Zhou_et_al_WXS, sheet = 1, startRow = 2, cols = 1:33)
# length(unique(zhou_et_al$Tumor_Sample_Barcode))

# Create a new column with Unique_Patient_Identifier for comparison with clinical data
zhou_et_al <- zhou_et_al %>%
  mutate(Unique_Patient_Identifier = sub("-.*", "", Tumor_Sample_Barcode))

zhou_et_al_maf <- preload_maf(maf = zhou_et_al, refset = "ces.refset.hg19", keep_extra_columns = "Hugo_Symbol")


### Clinical data
Zhou_et_al_source_data <- "zhou_et_al_source_data.xlsx"

if (!file.exists(Zhou_et_al_source_data)) {
  zhou_et_al_source_data_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-25787-x/MediaObjects/41467_2021_25787_MOESM10_ESM.xlsx"
  download.file(zhou_et_al_source_data_url, destfile = Zhou_et_al_source_data, mode = "wb")
}

zhou_et_al_clinical <- read.xlsx(Zhou_et_al_source_data, sheet = 9, startRow = 1, cols = 1:17)

# Remove unwanted columns
keep_cols <- c("patientID", "Age", "TNM","Smoking", "TreatmentAfterSurgery", "Status", "OS", "DFS")
zhou_et_al_clinical <- zhou_et_al_clinical[, keep_cols]

# Rename Columns
names(zhou_et_al_clinical) <- c("Unique_Patient_Identifier", 
                                 "age", 
                                "staging",
                                 "smoking_status",
                                 "previous.therapeutic.treatment.for.SCLC",
                                 "Status.(at.time.of.last.follow-up)", 
                                 "overall_survival.(months)",
                                 "progression-free_survival.(months)")


# Clean up data into appropriate format

zhou_et_al_clinical <- zhou_et_al_clinical %>%
  mutate(
    smoking_status = recode(smoking_status, "Smoker" = "smoker", "Nonsmoker" = "non-smoker"),
    `Status.(at.time.of.last.follow-up)` = recode(`Status.(at.time.of.last.follow-up)`, `0` = "alive", `1` = "dead"),
    previous.therapeutic.treatment.for.SCLC = recode(previous.therapeutic.treatment.for.SCLC, "Y" = TRUE, "N" = FALSE)
  )

zhou_et_al_clinical <- zhou_et_al_clinical %>%
  filter(Unique_Patient_Identifier %in% unique(zhou_et_al_maf$Unique_Patient_Identifier))


# # GENIE
# 
# GENIE_all_file_path_txt <- "GENIE_v17/data_mutations_extended.txt"
# GENIE_all_file_path_clinical <- "GENIE_v17/data_clinical_sample.txt"
# genomic_info_file_path <- "GENIE_v17/genomic_information.txt"
# 
# 
# 
# GENIE_all_txt <- read_tsv(GENIE_all_file_path_txt)
# GENIE_all_clinical <- read_tsv(GENIE_all_file_path_clinical, skip = 4)
# GENIE_SCLC_clinical <- GENIE_all_clinical[GENIE_all_clinical$ONCOTREE_CODE == "SCLC",]
# 
# keep_cols <- c("PATIENT_ID", "SAMPLE_ID", "AGE_AT_SEQ_REPORT","ONCOTREE_CODE", "SEQ_ASSAY_ID")
# 
# GENIE_SCLC_clinical <- GENIE_SCLC_clinical[, keep_cols]
# 
# names(GENIE_SCLC_clinical) <- c("Unique_Patient_Identifier", 
#                                          "Tumor_Sample_Barcode", 
#                                          "age",
#                                          "ONCOTREE_CODE", "SEQ_ASSAY_ID")
# GENIE_SCLC_txt <- GENIE_all_txt[GENIE_all_txt$Tumor_Sample_Barcode %in% GENIE_SCLC_clinical$Tumor_Sample_Barcode, ]
# 
# # Define coverage for targeted data set
# # Load genomic information data
# genie_genomic_info <- read_tsv(genomic_info_file_path, show_col_types = FALSE)
# 
# 
# # Convert genomic information into a GRanges object and ensure there are no outside sequence bounds
# # genie_genomic_info_GRanges <- makeGRangesFromDataFrame(genie_genomic_info, keep.extra.columns = TRUE, ignore.strand = TRUE, start.field = "Start_Position", end.field = "End_Position")
# 
# # # Did not find anything or work
# # view_17 <- genie_genomic_info_GRanges[genie_genomic_info_GRanges@seqnames == 17, ]
# # colSums(is.na(view_17@ranges))  # gives count of NAs in each column
# # 
# # which(end(genie_genomic_info_GRanges) > seqlengths(BSgenome.Hsapiens.UCSC.hg19)[as.character(seqnames(genie_genomic_info_GRanges))])
# # genie_genomic_info_GRanges.original <- genie_genomic_info_GRanges
# # genie_genomic_info_GRanges.trimmed <- trim(genie_genomic_info_GRanges)
# # which(genie_genomic_info_GRanges.original != genie_genomic_info_GRanges.trimmed)
#   
# 
# # Known hg19 chromosome lengths
# hg19_lengths <- read_tsv("https://hgdownload.soe.ucsc.edu/goldenPath/hg19/bigZips/hg19.chrom.sizes", n_max = 23)
# current_header <- colnames(hg19_lengths)
# hg19_lengths <- rbind(current_header, hg19_lengths)
# colnames(hg19_lengths) <- c("Chromosome", "Length")
# hg19_lengths$Chromosome <- sub("^chr", "", hg19_lengths$Chromosome)
# 
# # Filter all out-of-bound ranges for Chromosome 17
# genie_genomic_info <- left_join(genie_genomic_info, hg19_lengths, by = "Chromosome")
# # nrow(genie_genomic_info)
# genie_genomic_info <- genie_genomic_info %>%
#   filter(Chromosome != "17" | End_Position <= Length)
# genie_genomic_info <- genie_genomic_info[!(genie_genomic_info$SEQ_ASSAY_ID %in% c("MSK-IMPACT341", "MSK-IMPACT410", "MSK-IMPACT468")), ]
# 
# 
# # Create a mapping between samples and their sequencing assays
# # This links each SCLC sample to its corresponding sequencing panel
# sample_to_assay_mapping <- GENIE_SCLC_clinical %>%
#   select(Unique_Patient_Identifier, Tumor_Sample_Barcode, SEQ_ASSAY_ID)
# 
# # For patients with multiple samples/assays, keep only the one with the most comprehensive coverage
#   # Define the patients with multiple assays
# patients_with_multiple_assays <- sample_to_assay_mapping %>%
#   group_by(Unique_Patient_Identifier) %>%
#   summarise(n_assays = n_distinct(SEQ_ASSAY_ID), .groups = 'drop') %>%
#   filter(n_assays > 1) %>%
#   pull(Unique_Patient_Identifier)
# 
# for(patient_id in patients_with_multiple_assays) {
#   # Get all assays used for this patient
#   patient_assays <- sample_to_assay_mapping %>%
#     filter(Unique_Patient_Identifier == patient_id) %>%
#     pull(SEQ_ASSAY_ID) %>%
#     unique()
#   
#   # Find which assay covers the most genes
#   assay_gene_counts <- genie_genomic_info %>%
#     filter(SEQ_ASSAY_ID %in% patient_assays) %>%
#     distinct(SEQ_ASSAY_ID, Hugo_Symbol) %>%
#     count(SEQ_ASSAY_ID) %>%
#     arrange(desc(n))
#   
#   # Get the assay with the maximum number of genes
#   best_assay <- assay_gene_counts$SEQ_ASSAY_ID[1]
#   
#   # Remove samples using less comprehensive assays for this patient
#   sample_to_assay_mapping <- sample_to_assay_mapping %>%
#     filter(!(Unique_Patient_Identifier == patient_id & SEQ_ASSAY_ID != best_assay))
# }
# 
# # Add sequencing assay information to the SCLC MAF data
# GENIE_SCLC_txt_with_assay <- GENIE_SCLC_txt %>%
#   left_join(sample_to_assay_mapping, by = "Tumor_Sample_Barcode") %>%
#   # Remove any samples that don't have assay information
#   filter(!is.na(SEQ_ASSAY_ID))
# 
# # Preload the MAF data
# # This step prepares the data for cancer effects analysis
# GENIE_SCLC_maf_preloaded <- cancereffectsizeR::preload_maf(
#   maf = GENIE_SCLC_txt_with_assay,
#   refset = "ces.refset.hg19", 
#   sample_col = "Unique_Patient_Identifier", keep_extra_columns = c("Hugo_Symbol", "End_Position")
# )
# 
# 
# # Filter out germline variants and repetitive regions (except COSMIC-annotated ones)
# GENIE_SCLC_maf_preloaded <- GENIE_SCLC_maf_preloaded[
#   germline_variant_site == FALSE
# ][
#   repetitive_region == FALSE | cosmic_site_tier %in% 1:3
# ]
# 
# # Get list of unique assays used in our SCLC samples for subsequent loading into CESA object
# unique_assays <- unique(sample_to_assay_mapping$SEQ_ASSAY_ID)


#----------Import Data from Li et al., 2025 (Commentary Target)-----------------------

Li_et_al <- "li_et_al.xlsx"

if (!file.exists(Li_et_al)) {
  li_et_al_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41598-025-00222-z/MediaObjects/41598_2025_222_MOESM2_ESM.xlsx"
  download.file(li_et_al_url, destfile = Li_et_al, mode = "wb")
}

# Paper signature attribution biological weights

paper_signatures <- read.xlsx(Li_et_al, sheet = 3, startRow = 2, cols = 1:6)
paper_signatures$Sample <- sub(".*_", "", paper_signatures$Sample)
setnames(paper_signatures, old = "Sample", new = "Unique_Patient_Identifier", skip_absent=TRUE)
paper_signatures[!(paper_signatures$Unique_Patient_Identifier %in% george_et_al_clinical$Unique_Patient_Identifier), ]
george_et_al_clinical$Unique_Patient_Identifier
nrow(paper_signatures)

## WXS data
cox_regression_df <- read.xlsx(Li_et_al, sheet = 6, cols = 1:4, rows = 2:11963) %>%
  select(-Mutation_rate) %>%
  mutate(
    Univariate_Cox_fdr = p.adjust(Univariate_Cox_P, method = "BH"),
    Multivariate_Cox_fdr = p.adjust(Multivariate_Cox_P, method = "BH")
  )

top_x_genes <- 100 #Want to check the top 100 cancer effect genes for each ces_analysis

#-------- Import optimal cutoff 1000 bootstrap results -----
SBS13_TMB_optimal_cutoff_paper <- readRDS("SBS13_TMB_optimal_cutoff_paper.rds")
SBS13_TMB_optimal_cutoff_paper_mean <- round(SBS13_TMB_optimal_cutoff_paper[[2]], digits = 2)

SBS4_TMB_optimal_cutoff_paper <- readRDS("SBS4_TMB_optimal_cutoff_paper.rds")
SBS4_TMB_optimal_cutoff_paper_mean <- round(SBS4_TMB_optimal_cutoff_paper[[2]], 2)

SBS13_survival_optimal_cutoff_paper <- readRDS("SBS13_survival_optimal_cutoff_paper.rds")
SBS13_survival_optimal_cutoff_paper_mean <- round(SBS13_survival_optimal_cutoff_paper[[2]], 2)

SBS4_survival_optimal_cutoff_paper <- readRDS("SBS4_survival_optimal_cutoff_paper.rds")
SBS4_survival_optimal_cutoff_paper_mean <- round(SBS4_survival_optimal_cutoff_paper[[2]], 2)

# SBS13_TMB_optimal_cutoff_external <- readRDS("SBS13_TMB_optimal_cutoff_external.rds")
# SBS4_TMB_optimal_cutoff_external <- readRDS("SBS4_TMB_optimal_cutoff_external.rds")
# SBS13_survival_optimal_cutoff_external <- readRDS("SBS13_survival_optimal_cutoff_external.rds")
# SBS4_survival_optimal_cutoff_external <- readRDS("SBS4_survival_optimal_cutoff_external.rds")

