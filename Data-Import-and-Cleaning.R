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

