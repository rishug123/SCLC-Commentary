library(tidyverse)
library(data.table)
library(cancereffectsizeR)
library(maxstat)
library(survival)
library(ces.refset.hg19)

## -----------Import the cesa objects
cesa_paper <- load_cesa("cesa_paper_new.rds")
cesa_external <- load_cesa("cesa_external.rds")
cesa_all <- load_cesa("cesa_all.rds")

## -----------clear the assigned signatures
cesa_paper <- clear_trinuc_rates_and_signatures(cesa_paper)
cesa_external <- clear_trinuc_rates_and_signatures(cesa_external)
cesa_all <- clear_trinuc_rates_and_signatures(cesa_all)


## -----------extra code

#Define exclusions and assign mutational signatures
cosmic_confirmed_sig_for_sclc <- c(1, 4, 5, 15) # from https://cancer.sanger.ac.uk/signatures/signatures_v2/
reported_rel_sig_sclc <- c(1, 2, 3, 4, 5, 6, 13, 15, 16, 24, 29, 39, 40, 60, 92) # from https://www.nature.com/articles/s41586-024-07177-7 found to be present in n>5 cases and https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.891938/full highly confident signatures were identified with de-novo extraction.
#combining above lists for untreated
signature_untreated <- sort(unique(c(cosmic_confirmed_sig_for_sclc, reported_rel_sig_sclc)))
signature_untreated <- paste0("SBS", signature_untreated)
artifacts <- c("SBS27", "SBS43", "SBS45", "SBS46", "SBS47", "SBS48", "SBS49", "SBS50", "SBS51", 
               "SBS52", "SBS53", "SBS54", "SBS55", "SBS56", "SBS57", "SBS58", "SBS59", "SBS60", "SBS95")
#adding these for treated
treatment_signatures <- c("SBS11", "SBS31", "SBS32", "SBS35", "SBS86", "SBS87", "SBS90") # Note 
signature_treated <- c(signature_untreated, treatment_signatures)

sig_set <- c(ces.refset.hg19$signatures$COSMIC_v3.2$meta$Signature)
untreated_exclusions <- sig_set[!sig_set %in% c(signature_untreated, artifacts)]
treated_exclusions <- sig_set[!sig_set %in% c(signature_treated, artifacts)]

modified_signature_set <- ces.refset.hg19$signatures$COSMIC_v3.2
modified_signature_set$meta <- modified_signature_set$meta %>% 
  mutate(Likely_Artifact = ifelse(Signature == "SBS60", FALSE, Likely_Artifact))



filter_treated_untreated_with_maf <- function(samples_df, maf_df) {
  
  # Identify treated and untreated patient IDs
  untreated_ids <- samples_df %>%
    filter(`chemotherapy.(yes/no)` == "no" & previous.therapeutic.treatment.for.SCLC == FALSE) %>%
    pull(Unique_Patient_Identifier)
  
  treated_ids <- samples_df %>%
    filter(!(Unique_Patient_Identifier %in% untreated_ids)) %>%
    pull(Unique_Patient_Identifier)
  
  # Filter MAFs
  treated_maf <- maf_df %>% filter(Unique_Patient_Identifier %in% treated_ids)
  untreated_maf <- maf_df %>% filter(Unique_Patient_Identifier %in% untreated_ids)
  
  list(
    treated_ids = treated_ids,
    untreated_ids = untreated_ids,
    treated_maf = treated_maf,
    untreated_maf = untreated_maf
  )
}

## ----------bootstrap functions

bootstrap_analysis_SBS4_survival_optimal_cutoff_paper <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter clinical data
    clinical_data_paper <- cesa_paper$samples %>%
      filter(!is.na(`Status.(at.time.of.last.follow-up)`),
             !is.na(`overall_survival.(months)`)) %>%
      mutate(status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)) %>%
      select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)
    
    # Filter samples
    paper_filtered <- filter_treated_untreated_with_maf(cesa_paper$samples, cesa_paper$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = paper_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = paper_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_paper <- cesa_paper$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    # Merge clinical and signature data
    paper_optimal_cutoff <- biological_weights_paper %>%
      select(Unique_Patient_Identifier, SBS4, SBS13) %>%
      inner_join(clinical_data_paper, by = "Unique_Patient_Identifier")
    
    # Compute optimal cutoff using maxstat.test
    SBS4_paper_optimal_cutoff <- maxstat.test(
      Surv(`overall_survival.(months)`, status_numeric) ~ SBS4, data = paper_optimal_cutoff, smethod = "LogRank", 
      pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS4_paper_optimal_cutoff$estimate)
    
    cesa_paper <- clear_trinuc_rates_and_signatures(cesa_paper)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}
bootstrap_analysis_SBS13_survival_optimal_cutoff_paper <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter clinical data
    clinical_data_paper <- cesa_paper$samples %>%
      filter(!is.na(`Status.(at.time.of.last.follow-up)`),
             !is.na(`overall_survival.(months)`)) %>%
      mutate(status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)) %>%
      select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)
    
    # Filter samples
    paper_filtered <- filter_treated_untreated_with_maf(cesa_paper$samples, cesa_paper$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = paper_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = paper_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_paper <- cesa_paper$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    # Merge clinical and signature data
    paper_optimal_cutoff <- biological_weights_paper %>%
      select(Unique_Patient_Identifier, SBS4, SBS13) %>%
      inner_join(clinical_data_paper, by = "Unique_Patient_Identifier")
    
    # Compute optimal cutoff using maxstat.test
    SBS13_paper_optimal_cutoff <- maxstat.test(
      Surv(`overall_survival.(months)`, status_numeric) ~ SBS13, data = paper_optimal_cutoff, smethod = "LogRank", 
      pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS13_paper_optimal_cutoff$estimate)
    
    cesa_paper <- clear_trinuc_rates_and_signatures(cesa_paper)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}

bootstrap_analysis_SBS4_survival_optimal_cutoff_external <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter clinical data
    clinical_data_external <- cesa_external$samples %>%
      filter(!is.na(`Status.(at.time.of.last.follow-up)`),
             !is.na(`overall_survival.(months)`)) %>%
      mutate(status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)) %>%
      select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)
    
    # Filter samples
    external_filtered <- filter_treated_untreated_with_maf(cesa_external$samples, cesa_external$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_external <- trinuc_mutation_rates(
      cesa_external,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = external_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    # no samples seem to be untreated
    # cesa_external <- trinuc_mutation_rates(
    #   cesa_external,
    #   signature_extractor = "MutationalPatterns",
    #   signature_set = signature_set,
    #   signature_exclusions = untreated_exclusions,
    #   samples = external_filtered$untreated_ids,
    #   bootstrap_mutations = TRUE
    # )
    
    # Extract signature weights
    biological_weights_external <- cesa_external$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    # Merge clinical and signature data
    external_optimal_cutoff <- biological_weights_external %>%
      select(Unique_Patient_Identifier, SBS4, SBS13) %>%
      inner_join(clinical_data_external, by = "Unique_Patient_Identifier")
    
    # Compute optimal cutoff using maxstat.test
    SBS4_external_optimal_cutoff <- maxstat.test(
      Surv(`overall_survival.(months)`, status_numeric) ~ SBS4, data = external_optimal_cutoff, smethod = "LogRank", 
      pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS4_external_optimal_cutoff$estimate)
    
    cesa_external <- clear_trinuc_rates_and_signatures(cesa_external)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}
bootstrap_analysis_SBS13_survival_optimal_cutoff_external <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter clinical data
    clinical_data_external <- cesa_external$samples %>%
      filter(!is.na(`Status.(at.time.of.last.follow-up)`),
             !is.na(`overall_survival.(months)`)) %>%
      mutate(status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)) %>%
      select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)
    
    # Filter samples
    external_filtered <- filter_treated_untreated_with_maf(cesa_external$samples, cesa_external$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_external <- trinuc_mutation_rates(
      cesa_external,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = external_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    # all treated it seems like?
    # cesa_external <- trinuc_mutation_rates(
    #   cesa_external,
    #   signature_extractor = "MutationalPatterns",
    #   signature_set = signature_set,
    #   signature_exclusions = untreated_exclusions,
    #   samples = external_filtered$untreated_ids,
    #   bootstrap_mutations = TRUE
    # )
    
    # Extract signature weights
    biological_weights_external <- cesa_external$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    # Merge clinical and signature data
    external_optimal_cutoff <- biological_weights_external %>%
      select(Unique_Patient_Identifier, SBS4, SBS13) %>%
      inner_join(clinical_data_external, by = "Unique_Patient_Identifier")
    
    # Compute optimal cutoff using maxstat.test
    SBS13_external_optimal_cutoff <- maxstat.test(
      Surv(`overall_survival.(months)`, status_numeric) ~ SBS13, data = external_optimal_cutoff, smethod = "LogRank", 
      pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS13_external_optimal_cutoff$estimate)
    
    cesa_external <- clear_trinuc_rates_and_signatures(cesa_external)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}

bootstrap_analysis_SBS4_survival_optimal_cutoff_all <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter clinical data
    clinical_data_all <- cesa_all$samples %>%
      filter(!is.na(`Status.(at.time.of.last.follow-up)`),
             !is.na(`overall_survival.(months)`)) %>%
      mutate(status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)) %>%
      select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)
    
    # Filter samples
    all_filtered <- filter_treated_untreated_with_maf(cesa_all$samples, cesa_all$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = all_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = all_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_all <- cesa_all$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    # Merge clinical and signature data
    all_optimal_cutoff <- biological_weights_all %>%
      select(Unique_Patient_Identifier, SBS4, SBS13) %>%
      inner_join(clinical_data_all, by = "Unique_Patient_Identifier")
    
    # Compute optimal cutoff using maxstat.test
    SBS4_all_optimal_cutoff <- maxstat.test(
      Surv(`overall_survival.(months)`, status_numeric) ~ SBS4, data = all_optimal_cutoff, smethod = "LogRank", 
      pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS4_all_optimal_cutoff$estimate)
    
    cesa_all <- clear_trinuc_rates_and_signatures(cesa_all)
  }
  
  return(list(bs_samples, mean(bs_samples)))
}
bootstrap_analysis_SBS13_survival_optimal_cutoff_all <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter clinical data
    clinical_data_all <- cesa_all$samples %>%
      filter(!is.na(`Status.(at.time.of.last.follow-up)`),
             !is.na(`overall_survival.(months)`)) %>%
      mutate(status_numeric = ifelse(`Status.(at.time.of.last.follow-up)` == "dead", 1, 0)) %>%
      select(Unique_Patient_Identifier, status_numeric, `overall_survival.(months)`)
    
    # Filter samples
    all_filtered <- filter_treated_untreated_with_maf(cesa_all$samples, cesa_all$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = all_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = all_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_all <- cesa_all$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    # Merge clinical and signature data
    all_optimal_cutoff <- biological_weights_all %>%
      select(Unique_Patient_Identifier, SBS4, SBS13) %>%
      inner_join(clinical_data_all, by = "Unique_Patient_Identifier")
    
    # Compute optimal cutoff using maxstat.test
    SBS13_all_optimal_cutoff <- maxstat.test(
      Surv(`overall_survival.(months)`, status_numeric) ~ SBS13, data = all_optimal_cutoff, smethod = "LogRank", 
      pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS13_all_optimal_cutoff$estimate)
    
    cesa_all <- clear_trinuc_rates_and_signatures(cesa_all)
  }
  
  return(list(bs_samples, mean(bs_samples)))
}

bootstrap_analysis_SBS4_TMB_optimal_cutoff_paper <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter samples
    paper_filtered <- filter_treated_untreated_with_maf(cesa_paper$samples, cesa_paper$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = paper_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = paper_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_paper <- cesa_paper$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    paper_optimal_cutoff <- biological_weights_paper %>%
      select(Unique_Patient_Identifier, total_snvs, SBS4, SBS13) 
    
    # Compute optimal cutoff using maxstat.test
    SBS4_paper_optimal_cutoff <- maxstat.test(
      Surv(total_snvs) ~ SBS4, data = paper_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS4_paper_optimal_cutoff$estimate)
    
    cesa_paper <- clear_trinuc_rates_and_signatures(cesa_paper)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}
bootstrap_analysis_SBS13_TMB_optimal_cutoff_paper <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter samples
    paper_filtered <- filter_treated_untreated_with_maf(cesa_paper$samples, cesa_paper$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = paper_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_paper <- trinuc_mutation_rates(
      cesa_paper,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = paper_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_paper <- cesa_paper$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    paper_optimal_cutoff <- biological_weights_paper %>%
      select(Unique_Patient_Identifier, total_snvs, SBS4, SBS13) 
    
    # Compute optimal cutoff using maxstat.test
    SBS13_paper_optimal_cutoff <- maxstat.test(
      Surv(total_snvs) ~ SBS13, data = paper_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS13_paper_optimal_cutoff$estimate)
    
    cesa_paper <- clear_trinuc_rates_and_signatures(cesa_paper)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}

bootstrap_analysis_SBS4_TMB_optimal_cutoff_external <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter samples
    external_filtered <- filter_treated_untreated_with_maf(cesa_external$samples, cesa_external$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_external <- trinuc_mutation_rates(
      cesa_external,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = external_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    # cesa_external <- trinuc_mutation_rates(
    #   cesa_external,
    #   signature_extractor = "MutationalPatterns",
    #   signature_set = signature_set,
    #   signature_exclusions = untreated_exclusions,
    #   samples = external_filtered$untreated_ids,
    #   bootstrap_mutations = TRUE
    # )
    
    # Extract signature weights
    biological_weights_external <- cesa_external$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    external_optimal_cutoff <- biological_weights_external %>%
      select(Unique_Patient_Identifier, total_snvs, SBS4, SBS13) 
    
    # Compute optimal cutoff using maxstat.test
    SBS4_external_optimal_cutoff <- maxstat.test(
      Surv(total_snvs) ~ SBS4, data = external_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS4_external_optimal_cutoff$estimate)
    
    cesa_external <- clear_trinuc_rates_and_signatures(cesa_external)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}
bootstrap_analysis_SBS13_TMB_optimal_cutoff_external <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter samples
    external_filtered <- filter_treated_untreated_with_maf(cesa_external$samples, cesa_external$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_external <- trinuc_mutation_rates(
      cesa_external,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = external_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    # cesa_external <- trinuc_mutation_rates(
    #   cesa_external,
    #   signature_extractor = "MutationalPatterns",
    #   signature_set = signature_set,
    #   signature_exclusions = untreated_exclusions,
    #   samples = external_filtered$untreated_ids,
    #   bootstrap_mutations = TRUE
    # )
    
    # Extract signature weights
    biological_weights_external <- cesa_external$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    external_optimal_cutoff <- biological_weights_external %>%
      select(Unique_Patient_Identifier, total_snvs, SBS4, SBS13) 
    
    # Compute optimal cutoff using maxstat.test
    SBS13_external_optimal_cutoff <- maxstat.test(
      Surv(total_snvs) ~ SBS13, data = external_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS13_external_optimal_cutoff$estimate)
    
    cesa_external <- clear_trinuc_rates_and_signatures(cesa_external)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}

bootstrap_analysis_SBS4_TMB_optimal_cutoff_all <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter samples
    all_filtered <- filter_treated_untreated_with_maf(cesa_all$samples, cesa_all$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = all_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = all_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_all <- cesa_all$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    all_optimal_cutoff <- biological_weights_all %>%
      select(Unique_Patient_Identifier, total_snvs, SBS4, SBS13) 
    
    # Compute optimal cutoff using maxstat.test
    SBS4_all_optimal_cutoff <- maxstat.test(
      Surv(total_snvs) ~ SBS4, data = all_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS4_all_optimal_cutoff$estimate)
    
    cesa_all <- clear_trinuc_rates_and_signatures(cesa_all)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}
bootstrap_analysis_SBS13_TMB_optimal_cutoff_all <- function(num_bs) {
  bs_samples <- c()
  
  for (i in 1:num_bs) {
    
    # Filter samples
    all_filtered <- filter_treated_untreated_with_maf(cesa_all$samples, cesa_all$maf)
    signature_set <- modified_signature_set
    
    # Add trinucleotide mutation rates for treated and untreated
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = treated_exclusions,
      samples = all_filtered$treated_ids, 
      bootstrap_mutations = TRUE
    )
    
    cesa_all <- trinuc_mutation_rates(
      cesa_all,
      signature_extractor = "MutationalPatterns",
      signature_set = signature_set,
      signature_exclusions = untreated_exclusions,
      samples = all_filtered$untreated_ids,
      bootstrap_mutations = TRUE
    )
    
    # Extract signature weights
    biological_weights_all <- cesa_all$mutational_signatures$biological_weights[group_avg_blended == FALSE]
    
    all_optimal_cutoff <- biological_weights_all %>%
      select(Unique_Patient_Identifier, total_snvs, SBS4, SBS13) 
    
    # Compute optimal cutoff using maxstat.test
    SBS13_all_optimal_cutoff <- maxstat.test(
      Surv(total_snvs) ~ SBS13, data = all_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS13_all_optimal_cutoff$estimate)
    
    cesa_all <- clear_trinuc_rates_and_signatures(cesa_all)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}

## ---------testing to see if it works

# test1 <- bootstrap_analysis_SBS4_survival_optimal_cutoff_paper(5)
# test2 <- bootstrap_analysis_SBS4_TMB_optimal_cutoff_external(5)

## ---------what i would like to be processed with the cluster computer:

num_bs <- 1000 #Please set to whatever number is good

SBS4_survival_optimal_cutoff_paper <- bootstrap_analysis_SBS4_survival_optimal_cutoff_paper(num_bs)
SBS13_survival_optimal_cutoff_paper <- bootstrap_analysis_SBS13_survival_optimal_cutoff_paper(num_bs)
SBS4_TMB_optimal_cutoff_paper <- bootstrap_analysis_SBS4_TMB_optimal_cutoff_paper(num_bs)
SBS13_TMB_optimal_cutoff_paper <- bootstrap_analysis_SBS13_TMB_optimal_cutoff_paper(num_bs)

saveRDS(SBS4_survival_optimal_cutoff_paper, file = "SBS4_survival_optimal_cutoff_paper.rds")
saveRDS(SBS13_survival_optimal_cutoff_paper, file = "SBS13_survival_optimal_cutoff_paper.rds")
saveRDS(SBS4_TMB_optimal_cutoff_paper, file = "SBS4_TMB_optimal_cutoff_paper.rds")
saveRDS(SBS13_TMB_optimal_cutoff_paper, file = "SBS13_TMB_optimal_cutoff_paper.rds")

## ---------Extra for external and combined

# SBS4_survival_optimal_cutoff_external <- bootstrap_analysis_SBS4_survival_optimal_cutoff_external(num_bs)
# SBS13_survival_optimal_cutoff_external <- bootstrap_analysis_SBS13_survival_optimal_cutoff_external(num_bs)
# 
# SBS4_survival_optimal_cutoff_all <- bootstrap_analysis_SBS4_survival_optimal_cutoff_all(num_bs)
# SBS13_survival_optimal_cutoff_all <- bootstrap_analysis_SBS13_survival_optimal_cutoff_all(num_bs)
# 
# SBS4_TMB_optimal_cutoff_external <- bootstrap_analysis_SBS4_TMB_optimal_cutoff_external(num_bs)
# SBS13_TMB_optimal_cutoff_external <- bootstrap_analysis_SBS13_TMB_optimal_cutoff_external(num_bs)
# 
# SBS4_TMB_optimal_cutoff_all <- bootstrap_analysis_SBS4_TMB_optimal_cutoff_all(num_bs)
# SBS13_TMB_optimal_cutoff_all <- bootstrap_analysis_SBS13_TMB_optimal_cutoff_all(num_bs)


