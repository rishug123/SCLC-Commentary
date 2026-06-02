# Import external required packages
library(BiocManager)
library(cancereffectsizeR)
library(data.table)
library(tidyverse)
library(dplyr)
library(readr)
library(ggpmisc)
library(openxlsx)
library(readxl)
library(survminer)
library(survival)
library("maxstat")
library(BSgenome.Hsapiens.UCSC.hg19)
library(cowplot)
library(cutpointr)
library(gtools)
library(scales)
library(grid)
library(gridExtra)
library(contsurvplot)
library(pammtools)
library(riskRegression)
library(rms)
library(lmtest)
library(openxlsx2)


setwd("files")

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
keep_cols <- c("Sample-ID", "age", "sex", "ethnicity", "stage_UICC", "smoking_status", "previous.therapeutic.treatment.for.SCLC", "chemotherapy.(yes/no)","Status.(at.time.of.last.follow-up)", "progression-free_survival.(months)", "overall_survival.(months)")
george_et_al_clinical <- george_et_al_clinical[, keep_cols]

# Rename columns
setnames(george_et_al_clinical, old = "Sample-ID", new = "Unique_Patient_Identifier", skip_absent=TRUE)
setnames(george_et_al_clinical, old = "stage_UICC", new = "staging", skip_absent=TRUE)


george_et_al_clinical <- george_et_al_clinical[
  george_et_al_clinical$Unique_Patient_Identifier %in% unique(george_et_al_maf$Unique_Patient_Identifier), ]

george_et_al_clinical <- george_et_al_clinical %>%
  mutate(
    previous.therapeutic.treatment.for.SCLC = recode(previous.therapeutic.treatment.for.SCLC, "untreated"= FALSE, "relapse" = TRUE),
    smoking_status = if_else(is.na(smoking_status), "non-smoker", "smoker"),
    staging = toupper(trimws(staging)),
    ethnicity = tolower(trimws(ethnicity))
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
keep_cols <- c("Unique_Patient_Identifier", "Age", "Gender", "Race", "UICC.stage", "Smoke", "Months", "Survival.Status", "chemotherapy-exposed.biopsy")
jiang_et_al_clinical <- jiang_et_al_clinical[, keep_cols]

# Rename Columns
names(jiang_et_al_clinical) <- c("Unique_Patient_Identifier",
                                 "age",
                                 "sex", "ethnicity", "staging",
                                 "smoking_status",
                                 "overall_survival.(months)",
                                 "Status.(at.time.of.last.follow-up)", "chemotherapy.(yes/no)")

# Clean up data into appropriate format

jiang_et_al_clinical <- jiang_et_al_clinical %>%
  mutate(
    sex = recode(sex, `1` = "male", `2` = "female"),
    smoking_status = recode(smoking_status, `1` = "smoker", `2` = "non-smoker"),
    `Status.(at.time.of.last.follow-up)` = recode(`Status.(at.time.of.last.follow-up)`, `0` = "dead", `1` = "alive"),
    `chemotherapy.(yes/no)` = recode(`chemotherapy.(yes/no)`, `Chemo` = "yes", `Naive` = "no"),
    staging = if_else(grepl("/", staging, fixed = TRUE), NA_character_, toupper(trimws(staging))),
    ethnicity = tolower(trimws(ethnicity))
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
    previous.therapeutic.treatment.for.SCLC = recode(previous.therapeutic.treatment.for.SCLC, "Y" = TRUE, "N" = FALSE),
    staging = toupper(trimws(staging)),
    ethnicity = "chinese"
  )

zhou_et_al_clinical <- zhou_et_al_clinical %>%
  filter(Unique_Patient_Identifier %in% unique(zhou_et_al_maf$Unique_Patient_Identifier))


## Song et al., 2021 (Cell Death & Disease)

### WXS
Song_et_al <- "song_et_al_2021.xlsx"

if (!file.exists(Song_et_al)) {
  song_et_al_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41419-021-03754-0/MediaObjects/41419_2021_3754_MOESM2_ESM.xlsx"
  download.file(song_et_al_url, destfile = Song_et_al, mode = "wb")
}

song_et_al <- as.data.frame(read_excel(Song_et_al, sheet = "Table S3"))

setnames(song_et_al, old = "Patient_ID", new = "Unique_Patient_Identifier", skip_absent = TRUE)
setnames(song_et_al, old = "Tumor_Seq_Allele1", new = "Tumor_Seq_Allele", skip_absent = TRUE)

song_et_al <- song_et_al %>%
  filter(!is.na(Hugo_Symbol), !is.na(Start_Position)) %>%
  mutate(Chromosome = paste0("chr", Chromosome))

song_et_al_maf <- preload_maf(maf = song_et_al, refset = "ces.refset.hg19", tumor_allele_col = "Tumor_Seq_Allele", keep_extra_columns = maf_format)


### Clinical data
song_et_al_clinical <- as.data.frame(read_excel(Song_et_al, sheet = "Table S1"))
song_et_al_clinical <- song_et_al_clinical[song_et_al_clinical$Cancer == "SCLC", ]

keep_cols <- c("Patient_ID", "Smoking", "Stage", "OS (months)", "Vital status")
song_et_al_clinical <- song_et_al_clinical[, keep_cols]

names(song_et_al_clinical) <- c("Unique_Patient_Identifier",
                                "smoking_status",
                                "staging",
                                "overall_survival.(months)",
                                "Status.(at.time.of.last.follow-up)")

song_et_al_clinical <- song_et_al_clinical %>%
  dplyr::mutate(smoking_status = dplyr::case_when(
    smoking_status == 1 ~ "smoker",
    smoking_status == 0 ~ "non-smoker",
    TRUE ~ NA_character_
  ),
  `Status.(at.time.of.last.follow-up)` = dplyr::case_when(
    `Status.(at.time.of.last.follow-up)` == 1 ~ "dead",
    `Status.(at.time.of.last.follow-up)` == 0 ~ "alive",
    TRUE ~ NA_character_
  ),
  staging = dplyr::case_when(
    is.na(staging) | staging == "NA" ~ NA_character_,
    TRUE ~ toupper(trimws(staging))
  ),
  ethnicity = "chinese") %>%
  filter(Unique_Patient_Identifier %in% unique(song_et_al_maf$Unique_Patient_Identifier))


## Chen et al., 2021 (Nature Communications) - multi-region, keep T1 only

### WXS
Chen_et_al_WXS <- "chen_et_al_2021_maf.xlsx"

if (!file.exists(Chen_et_al_WXS)) {
  chen_et_al_WXS_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-26821-8/MediaObjects/41467_2021_26821_MOESM4_ESM.xlsx"
  download.file(chen_et_al_WXS_url, destfile = Chen_et_al_WXS, mode = "wb")
}

# Barcodes are "P01-T1-T", "P01-T2-T", ... - keep T1 only to avoid TMB inflation
chen_et_al <- as.data.frame(read_excel(Chen_et_al_WXS, sheet = "data 2")) %>%
  filter(grepl("-T1-", tumor_name)) %>%
  mutate(
    Hugo_Symbol            = gene,
    Chromosome             = paste0("chr", chr),
    Start_Position         = as.integer(start),
    End_Position           = as.integer(end),
    Reference_Allele       = ref_allele,
    Tumor_Seq_Allele       = alt_allele,
    Variant_Classification = if_else(!is.na(exonicfunc.knowngene) & exonicfunc.knowngene != ".",
                                     exonicfunc.knowngene, func.knowngene),
    Unique_Patient_Identifier = paste0("chen_", sub("-.*", "", tumor_name)),
    Tumor_Sample_Barcode      = Unique_Patient_Identifier
  ) %>%
  filter(!is.na(Hugo_Symbol), !is.na(Start_Position))

chen_et_al_maf <- preload_maf(maf = chen_et_al, refset = "ces.refset.hg19", tumor_allele_col = "Tumor_Seq_Allele", keep_extra_columns = maf_format)


### Clinical data
Chen_et_al_clinical <- "chen_et_al_2021_clinical.xlsx"

if (!file.exists(Chen_et_al_clinical)) {
  chen_et_al_clinical_url <- "https://static-content.springer.com/esm/art%3A10.1038%2Fs41467-021-26821-8/MediaObjects/41467_2021_26821_MOESM3_ESM.xlsx"
  download.file(chen_et_al_clinical_url, destfile = Chen_et_al_clinical, mode = "wb")
}

chen_et_al_clinical <- as.data.frame(read_excel(Chen_et_al_clinical, sheet = "data 1", skip = 2)) %>%
  filter(grepl("^P[0-9]+$", `New ID`))

keep_cols <- c("New ID", "Smoking", "TNM Stage", "OS (months)", "Vital status")
chen_et_al_clinical <- chen_et_al_clinical[, keep_cols]

names(chen_et_al_clinical) <- c("Unique_Patient_Identifier",
                                "smoking_status",
                                "staging",
                                "overall_survival.(months)",
                                "Status.(at.time.of.last.follow-up)")

chen_et_al_clinical <- chen_et_al_clinical  %>%
  dplyr::mutate(
    smoking_status = dplyr::case_when(
      smoking_status == "Yes" ~ "smoker",
      smoking_status == "No" ~ "non-smoker",
      TRUE ~ NA_character_
    ),
    `Status.(at.time.of.last.follow-up)` = dplyr::case_when(
      `Status.(at.time.of.last.follow-up)` == "Died" ~ "dead",
      `Status.(at.time.of.last.follow-up)` == "Alive" ~ "alive",
      TRUE ~ NA_character_
    ),
    staging = toupper(trimws(staging)),
    ethnicity = NA_character_
  ) %>%
  mutate(
    Unique_Patient_Identifier = paste0("chen_", Unique_Patient_Identifier)
  ) %>%
  filter(Unique_Patient_Identifier %in% unique(chen_et_al_maf$Unique_Patient_Identifier))


## Wang et al., 2022 (Cancer Science)

### WXS
Wang_et_al_WXS <- "wang_et_al_2022_tableS4.xls"

if (!file.exists(Wang_et_al_WXS)) {
  wang_et_al_WXS_url <- "https://onlinelibrary.wiley.com/action/downloadSupplement?doi=10.1111%2Fcas.15606&file=cas15606-sup-0016-TableS4.xls"
  download.file(wang_et_al_WXS_url, destfile = Wang_et_al_WXS, mode = "wb")
}

wang_vc_map <- c(
  "nonsynonymous SNV"          = "Missense_Mutation",
  "stopgain"                   = "Nonsense_Mutation",
  "stoploss"                   = "Nonstop_Mutation",
  "startloss"                  = "Translation_Start_Site",
  "frameshift deletion"        = "Frame_Shift_Del",
  "frameshift insertion"       = "Frame_Shift_Ins",
  "frameshift substitution"    = "Frame_Shift_Sub",
  "nonframeshift deletion"     = "In_Frame_Del",
  "nonframeshift insertion"    = "In_Frame_Ins",
  "nonframeshift substitution" = "In_Frame_Del"
)

wang_et_al <- as.data.frame(read_excel(Wang_et_al_WXS, sheet = "Table S4A-SNVs-Chinese", skip = 2)) %>%
  mutate(
    Variant_Classification = if_else(ExonicFunc.refGene != "." & !is.na(ExonicFunc.refGene),
                                     wang_vc_map[ExonicFunc.refGene], "Splice_Site"),
    Start_Position = as.integer(Start),
    End_Position   = as.integer(End)
  ) %>%
  filter(!is.na(Gene.refGene), !is.na(Start_Position))

setnames(wang_et_al, old = "Sample ID", new = "Unique_Patient_Identifier", skip_absent = TRUE)
setnames(wang_et_al, old = "Chr",       new = "Chromosome",                skip_absent = TRUE)
setnames(wang_et_al, old = "Ref",       new = "Reference_Allele",          skip_absent = TRUE)
setnames(wang_et_al, old = "Alt",       new = "Tumor_Seq_Allele",          skip_absent = TRUE)
setnames(wang_et_al, old = "Gene.refGene", new = "Hugo_Symbol",            skip_absent = TRUE)

wang_et_al_maf <- preload_maf(maf = wang_et_al, refset = "ces.refset.hg19", tumor_allele_col = "Tumor_Seq_Allele", keep_extra_columns = maf_format)

## Clinical Data Unavailable

## Liu et al., 2024 (Cell)

if (!file.exists("liu_et_al_wxs_wb.rds") || !file.exists("liu_et_al_clinical_wb.rds")) {
  liu_et_al <- wb_load("liu_et_al_2024_tableS1.xlsb")
  
  liu_et_al_wxs_df      <- wb_to_df(liu_et_al, sheet = 3)
  liu_et_al_clinical_df <- wb_to_df(liu_et_al, sheet = 2)
  
  saveRDS(liu_et_al_wxs_df, "liu_et_al_wxs_wb.rds")
  saveRDS(liu_et_al_clinical_df, "liu_et_al_clinical_wb.rds")
  
} else {
  liu_et_al_wxs_df      <- readRDS("liu_et_al_wxs_wb.rds")
  liu_et_al_clinical_df <- readRDS("liu_et_al_clinical_wb.rds")
}

### WXS
liu_et_al <- liu_et_al_wxs_df %>%
  mutate(
    Unique_Patient_Identifier = paste0("liu_", sub("^T", "", Tumor_Sample_Barcode)),
    Tumor_Sample_Barcode      = Unique_Patient_Identifier
  ) %>%
  filter(!is.na(Hugo_Symbol), !is.na(Start_Position)) %>% dplyr::select(!Genome_Change)

liu_et_al_maf <- preload_maf(maf = liu_et_al, refset = "ces.refset.hg19",
                              chain_file = "hg38ToHg19.over.chain",
                              tumor_allele_col = "Tumor_Seq_Allele2",
                              keep_extra_columns = maf_format)

### Clinical data
keep_cols <- c("Sample.ID", "Age", "Gender", "TNM.Stage", "Smoking",
               "Status.", "Survial.(months)")
liu_et_al_clinical <- liu_et_al_clinical_df[, keep_cols]

names(liu_et_al_clinical) <- c("Unique_Patient_Identifier",
                                "age",
                                "sex",
                                "staging",
                                "smoking_status",
                                "Status.(at.time.of.last.follow-up)",
                                "overall_survival.(months)")

liu_et_al_clinical <- liu_et_al_clinical %>%
  mutate(
    Unique_Patient_Identifier = paste0("liu_", Unique_Patient_Identifier),
    sex = dplyr::case_when(
      sex == "M" ~ "male",
      sex == "F" ~ "female",
      TRUE ~ NA_character_
    ),
    smoking_status = dplyr::case_when(
      smoking_status == "yes"  ~ "smoker",
      smoking_status == "no"   ~ "non-smoker",
      TRUE ~ NA_character_
    ),
    `Status.(at.time.of.last.follow-up)` = tolower(trimws(`Status.(at.time.of.last.follow-up)`)),
    `overall_survival.(months)` = as.numeric(
      sub("^([0-9.]+).*", "\\1", as.character(`overall_survival.(months)`))
    ),
    staging = toupper(trimws(gsub("[0-9]", "", staging))),
    ethnicity = "chinese"
  ) %>%
  filter(Unique_Patient_Identifier %in% unique(liu_et_al_maf$Unique_Patient_Identifier))


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

# #-------- Import optimal cutoff 1000 bootstrap results -----
# SBS13_TMB_optimal_cutoff_paper <- readRDS("SBS13_TMB_optimal_cutoff_paper.rds")
# SBS13_TMB_optimal_cutoff_paper_mean <- round(SBS13_TMB_optimal_cutoff_paper[[2]], digits = 2)
# 
# SBS4_TMB_optimal_cutoff_paper <- readRDS("SBS4_TMB_optimal_cutoff_paper.rds")
# SBS4_TMB_optimal_cutoff_paper_mean <- round(SBS4_TMB_optimal_cutoff_paper[[2]], 2)
# 
# SBS13_survival_optimal_cutoff_paper <- readRDS("SBS13_survival_optimal_cutoff_paper.rds")
# SBS13_survival_optimal_cutoff_paper_mean <- round(SBS13_survival_optimal_cutoff_paper[[2]], 2)
# 
# SBS4_survival_optimal_cutoff_paper <- readRDS("SBS4_survival_optimal_cutoff_paper.rds")
# SBS4_survival_optimal_cutoff_paper_mean <- round(SBS4_survival_optimal_cutoff_paper[[2]], 2)

