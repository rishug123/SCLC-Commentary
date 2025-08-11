# Create a summary of all genomic and clinical data

summarize_maf_sources <- function(cesa_all) {
  # Extract unique sources
  sources <- unique(cesa_all$samples$maf_source)
  
  # Initialize results data.frame
  results <- data.frame(
    maf_source = character(),
    num_patients = integer(),
    num_with_survival = integer(),
    num_maf_rows = integer(),
    stringsAsFactors = FALSE
  )
  
  # Loop through each source
  for (source in sources) {
    pts <- cesa_all$samples %>% 
      filter(maf_source == source)
    
    # Count patients
    num_patients <- nrow(pts)
    
    # Count patients with survival
    num_with_survival <- sum(!is.na(pts$`overall_survival.(months)`))
    
    # MAF rows for these patients
    maf_rows <- cesa_all$maf %>%
      filter(Unique_Patient_Identifier %in% pts$Unique_Patient_Identifier)
    num_maf_rows <- nrow(maf_rows)
    
    # Append to results
    results <- rbind(
      results,
      data.frame(
        maf_source = source,
        num_patients = num_patients,
        num_with_survival = num_with_survival,
        num_maf_rows = num_maf_rows,
        stringsAsFactors = FALSE
      )
    )
  }
  
  return(results)
}

# ------------------Signature analysis-------------------

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

# Function to assign mutational signatures

assign_mutational_signatures <- function(signature_extractor_type) {
  if (!signature_extractor_type %in% c("MutationalPatterns", "deconstructSigs")) {
    stop("Invalid signature extractor type. Must be 'MutationalPatterns' or 'deconstructSigs'.")
  }
  
  paper_filtered <- filter_treated_untreated_with_maf(
    cesa_paper$samples, 
    cesa_paper$maf
    )
  
  external_filtered <- filter_treated_untreated_with_maf(
    cesa_external$samples,
    cesa_external$maf
    )
  
  all_filtered <- filter_treated_untreated_with_maf(
    cesa_all$samples, 
    cesa_all$maf
    )
  
  # Redefine SBS60 signature to Unknown from Potential sequencing artifact to include it in analysis
  manipulated_sig_set <- ces.refset.hg19$signatures$COSMIC_v3.2
  manipulated_sig_set$meta <- manipulated_sig_set$meta %>%
    mutate(
      Etiology = ifelse(Signature == "SBS60", "Unknown", Etiology),
      Likely_Artifact = ifelse(Signature == "SBS60", FALSE, Likely_Artifact)
    )
  artifact_names <- manipulated_sig_set$meta[Likely_Artifact == TRUE, Signature]
  signature_set <- manipulated_sig_set
  
  #Define exclusions and assign mutational signatures
  all_sigs <- rownames(ces.refset.hg19$signatures$COSMIC_v3.2$signatures)
  cosmic_confirmed_sig_for_sclc <- c(1, 4, 5, 15) # from https://cancer.sanger.ac.uk/signatures/signatures_v2/
  reported_rel_sig_sclc <- c(1, 2, 3, 4, 5, 6, 13, 15, 16, 24, 29, 39, 40, 60, 92) # from https://www.nature.com/articles/s41586-024-07177-7 found to be present in n>5 cases and https://www.frontiersin.org/journals/oncology/articles/10.3389/fonc.2022.891938/full highly confident signatures were identified with de-novo extraction.
  
  sclc_sigs <- sort(unique(c(cosmic_confirmed_sig_for_sclc, reported_rel_sig_sclc)))
  signature_untreated <- paste0("SBS", sclc_sigs)
  
  treatment_signatures <- c("SBS11", "SBS31", "SBS32", "SBS35", "SBS86", "SBS87", "SBS90") # Note 
  signature_treated <- c(signature_untreated, treatment_signatures)
  
  signature_exclusions_untreated <- setdiff(all_sigs, signature_untreated)
  signature_exclusions_treated <- setdiff(all_sigs, signature_treated)
  signature_exclusions_untreated <- setdiff(signature_exclusions_untreated, artifact_names)
  signature_exclusions_treated <- setdiff(signature_exclusions_treated, artifact_names)
  
  if (length(paper_filtered$treated_ids) > 0) {
    cesa_paper <<- trinuc_mutation_rates(
      cesa_paper, 
      signature_extractor = signature_extractor_type,
      signature_set = signature_set,
      signature_exclusions = signature_exclusions_treated,
      samples = paper_filtered$treated_ids
      )
  }
  if (length(paper_filtered$untreated_ids) > 0) {
    cesa_paper <<- trinuc_mutation_rates(
      cesa_paper, 
      signature_extractor = signature_extractor_type,
      signature_set = signature_set,
      signature_exclusions = signature_exclusions_untreated,
      samples = paper_filtered$untreated_ids
      )
  }
  
  if (length(external_filtered$treated_ids) > 0) {
    cesa_external <<- trinuc_mutation_rates(
      cesa_external,
      signature_extractor = signature_extractor_type,
      signature_set = signature_set,
      signature_exclusions = signature_exclusions_treated,
      samples = external_filtered$treated_ids
      )
  }
  if (length(external_filtered$untreated_ids) > 0) {
    cesa_external <<- trinuc_mutation_rates(
      cesa_external, 
      signature_extractor = signature_extractor_type,
      signature_set = signature_set,
      signature_exclusions = signature_exclusions_untreated,
      samples = external_filtered$untreated_ids
      )
  }
  
  if (length(all_filtered$treated_ids) > 0) {
    cesa_all <<- trinuc_mutation_rates(
      cesa_all, 
      signature_extractor = signature_extractor_type,
      signature_set = signature_set,
      signature_exclusions = signature_exclusions_treated,
      samples = all_filtered$treated_ids
      )
  }
  if (length(all_filtered$untreated_ids) > 0) {
    cesa_all <<- trinuc_mutation_rates(
      cesa_all, 
      signature_extractor = signature_extractor_type,
      signature_set = signature_set,
      signature_exclusions = signature_exclusions_untreated,
      samples = all_filtered$untreated_ids
      )
  }

}

# Determine the optimal binary for signature presence

bootstrap_analysis_SBS24_survival_optimal_cutoff_paper <- function(num_bs) {
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
      select(Unique_Patient_Identifier, SBS24) %>%
      inner_join(clinical_data_paper, by = "Unique_Patient_Identifier")
    
    # Compute optimal cutoff using maxstat.test
    SBS24_paper_optimal_cutoff <- maxstat.test(
      Surv(`overall_survival.(months)`, status_numeric) ~ SBS24, data = paper_optimal_cutoff, smethod = "LogRank", 
      pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS24_paper_optimal_cutoff$estimate)
    
    cesa_paper <- clear_trinuc_rates_and_signatures(cesa_paper)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}
bootstrap_analysis_SBS24_TMB_optimal_cutoff_paper <- function(num_bs) {
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
      select(Unique_Patient_Identifier, total_snvs, SBS24) 
    
    # Compute optimal cutoff using maxstat.test
    SBS24_paper_optimal_cutoff <- maxstat.test(
      Surv(total_snvs) ~ SBS24, data = paper_optimal_cutoff, smethod = "LogRank", pmethod = "HL", minprop = 0.01
    )
    
    # Append to bootstrap samples
    bs_samples <- c(bs_samples, SBS24_paper_optimal_cutoff$estimate)
    
    cesa_paper <- clear_trinuc_rates_and_signatures(cesa_paper)
  }
  
  # Compute mean of bootstrap estimates
  return(list(bs_samples, mean(bs_samples)))
}


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

plot_TMB_binary_boxplot <- function(data, signature_col, cutoff, label) {
  font <- 12
  num_digits_after_decimal <- 2
  
  # Step 1: Create 2-line x-axis labels
  data <- data %>%
    mutate(
      x_group = ifelse(.data[[signature_col]] >= cutoff,
                       paste0("High (\u2265 ", cutoff, ")"),
                       paste0("Low (< ", cutoff, ")"))
    )
  
  # Step 2: Count group sizes
  group_counts <- table(data$x_group)
  high_label <- paste0("High (\u2265 ", cutoff, ")")
  low_label  <- paste0("Low (< ", cutoff, ")")
  
  high_count <- ifelse(high_label %in% names(group_counts), group_counts[high_label], 0)
  low_count  <- ifelse(low_label %in% names(group_counts), group_counts[low_label], 0)
  
  # Calculate Wilcoxon p-value
  p_value <- wilcox.test(total_snvs ~ x_group, data = data)$p.value
  if (p_value > 0.01){
    p_label <- paste0("italic(P) == ", signif(p_value, num_digits_after_decimal))
  } else{
    p_value <- gsub("e(-?\\d+)", " %*% 10^{\\1}", formatC(p_value, format = "e", digits = num_digits_after_decimal))
    p_label <- paste0("italic(P) == ", p_value)
  }
  
  # Step 3: Color group labels (keep on one line)
  data <- data %>%
    mutate(
      color_group = case_when(
        x_group == high_label ~ paste0("High (n=", high_count, ")"),
        x_group == low_label  ~ paste0("Low (n=", low_count, ")")
      ),
      x_group = factor(x_group, levels = c(low_label, high_label)),
      color_group = factor(color_group, levels = c(
        paste0("Low (n=", low_count, ")"),
        paste0("High (n=", high_count, ")")
      ))
    )
  
  
  # Step 4: Create the base plot with corrected legend title
  base_plot <- ggboxplot(data,
                         x = "x_group",
                         y = "total_snvs",
                         color = "color_group",
                         palette = c("red", "blue"),
                         add = "jitter") +
    
    scale_y_continuous(
      trans = "log10",
      limits = c(10^0, 10^3),
      breaks = 10^seq(1, 5, by = 1),
      labels = function(x) parse(text = paste0("10^", log10(x)))
    ) +
    
    ylab("TMB") +
    xlab(paste0(signature_col, " signature status")) +
    
    # stat_compare_means(
    #   label = "p.format", size = font / 2.8, method = "wilcox.test", label.x = 0, label.y = 10^1
    # ) +
    annotate("text", x = 2, y = 10^0.75, label = p_label, parse = TRUE, size = font / 2.845) +
    
    theme_bw() +
    theme(
      plot.title = element_text(size = font, face = "bold", color = "black"),
      axis.title = element_text(size = font, face = "plain", color = "black"),
      axis.text = element_text(size = font),
      legend.title = element_text(size = font, face = "plain", color = "black"),  # Fixed: proper styling
      legend.text = element_text(size = font, face = "plain", color = "black"),
      legend.position = "top",
      plot.margin = margin(t = 1, r = 1, b = 1, l = 5)
    ) +
    
    # Add the legend title using labs()
    labs(color = paste0(signature_col, " status")) +
    
    # Format the legend with proper spacing
    guides(color = guide_legend(title = paste0(signature_col, " status\n"),
                                title.position = "top",
                                title.hjust = 0.5))
  
  return(base_plot)
}

plot_TMB_linear_regression <- function(biological_weights_table, dataset, signature, cutoff_SBS4 = SBS4_survival_optimal_cutoff_paper_mean, cutoff_SBS13 = SBS13_survival_optimal_cutoff_paper_mean) {
  
  font <- 12 # Master font size
  num_digits_after_decimal <- 2
  
  # Helper: Prepare table for each signature
  make_table <- function(sig) {
    data.table(
      Signature = paste0(sig, " (n=", nrow(biological_weights_table), ")"),
      Proportion_Signature = biological_weights_table[[sig]],
      log10_TMB = biological_weights_table$total_snvs
    )
  }
  
  table <- make_table(signature)
  
  # Plot builder
  build_plot <- function(data, sig_name, color, cutoff = NA, xmax = 0.8){
    
    model <- lm(log10_TMB ~ Proportion_Signature, data = data)
    coefs <- coef(model)
    summary_lm <- summary(model)
    r2_val <- signif(summary_lm$r.squared, num_digits_after_decimal)
    p_val  <- summary_lm$coefficients[2, 4]
    
    # Equation string
    eqn_string <- paste0(
      "italic(y) == ", signif(coefs[1], num_digits_after_decimal),
      if (coefs[2] >= 0) " + " else " - ",
      abs(signif(coefs[2], num_digits_after_decimal)), " * italic(x)"
    )
    
    # R² string
    r2_string <- paste0("italic(R)^2 == ", r2_val)
    
    # P-value string
    if (p_val > 0.01){
      p_label <- paste0("italic(P) == ", signif(p_val, num_digits_after_decimal))
    } else{
      p_value <- gsub("e(-?\\d+)", " %*% 10^{\\1}", formatC(p_val, format = "e", digits = num_digits_after_decimal))
      p_label <- paste0("italic(P) == ", p_value)
    }
    
    ggplot(data, aes(x = Proportion_Signature, y = log10_TMB)) +
      geom_point(size = 1.5, color = color) +
      geom_smooth(method = "lm", se = FALSE, color = color, size = 0.5) +
      scale_y_continuous(
        trans = "log10",
        limits = c(10^0, 10^3),
        breaks = 10^(1:3),
        labels = trans_format("log10", math_format(10^.x))
      ) +
      xlim(0, xmax) +
      labs(x = paste0(sig_name, " signature activity"), y = "TMB") +
      annotate("text", x = xmax, y = 10^1, label = eqn_string, parse = TRUE,
               hjust = 1.05, vjust = 0, size = font / 4, color = color) +
      annotate("text", x = xmax, y = 10^0.6, label = r2_string, parse = TRUE,
               hjust = 1.05, vjust = 0, size = font / 4, color = color) +
      annotate("text", x = xmax, y = 10^0.2, label = p_label, parse = TRUE,
               hjust = 1.05, vjust = 0, size = font / 4, color = color) +
      theme_bw(base_size = font) +
      theme(
        # panel.grid.minor = element_blank(),
        # panel.grid.major = element_blank(),
        axis.text = element_text(size = font * 0.8),
        axis.title = element_text(size = font),
        plot.title = element_blank()
      ) #   + { if (!is.na(cutoff)) geom_vline(xintercept = cutoff, linetype = "dashed", color = color) else NULL } 
    }
  
  # Pick xmax dynamically
  xmax <- if (signature == "SBS13") 0.4 else 0.8
  
  # Call builder
  color_map <- list(
    SBS4  = "#a6761d",
    SBS13 = "#7570b3",
    SBS5  = "#999999",
    SBS87 = "#1b9e77",
    SBS3  = "#e7298a",
    SBS24 = "#66a61e"
  )
  
  cutoff <- switch(signature,
                   SBS4 = cutoff_SBS4,
                   SBS13 = cutoff_SBS13,
                   NA)
  
  color <- color_map[[signature]]
  
  plot <- build_plot(table, signature, color, cutoff = cutoff, xmax = xmax)
  return(plot)
}


analyze_signature <- function(signature_name, dataset, biological_weight_table, survival_data, cutoff = NULL, SBS4_cutoff = SBS4_survival_optimal_cutoff_paper_mean, SBS13_cutoff = SBS13_survival_optimal_cutoff_paper_mean) {
  
  # Define styling constants
  font <- 12
  font_legend <- 9
  risk_table_number_size <- font/2.5
  
  # Helper function to create binary groups
  create_binary_groups <- function(threshold, threshold_label) {
    # More efficient filtering using merge instead of nested subsetting
    high_patients <- biological_weight_table[biological_weight_table[[signature_name]] >= threshold, "Unique_Patient_Identifier", drop = FALSE]
    low_patients <- biological_weight_table[biological_weight_table[[signature_name]] < threshold, "Unique_Patient_Identifier", drop = FALSE]
    
    # Use merge for more efficient joins
    high_group <- merge(survival_data, high_patients, by = "Unique_Patient_Identifier")
    low_group <- merge(survival_data, low_patients, by = "Unique_Patient_Identifier")
    
    # Add group labels
    high_group$group <- paste0("High (\u2265 ", threshold_label, ")")
    low_group$group <- paste0("Low (< ", threshold_label, ")")
    
    # Combine and set factor levels
    combined_df <- rbind(high_group, low_group)
    combined_df$group <- factor(combined_df$group, 
                                levels = c(paste0("Low (< ", threshold_label, ")"), 
                                           paste0("High (\u2265 ", threshold_label, ")")))
    
    return(list(
      combined_df = combined_df,
      high_group = high_group,
      low_group = low_group
    ))
  }
  
  # Helper function to create survival plot
  create_survival_plot <- function(surv_fit, surv_diff, data, threshold_label, dataset) {
    p_label <- paste0("italic(P) == ", signif(surv_diff$pvalue, 2))
    plot <- ggsurvplot(
      fit = surv_fit,
      data = data,
      xlim = c(0, 80),
      break.x.by = 20,
      ylab = "Overall survival",
      xlab = "Time (months)",
      pval = FALSE,
      risk.table = TRUE,
      risk.table.height = 0.20,
      legend.labs = c(paste0(signature_name, " low (< ", threshold_label, ")"), 
                      paste0(signature_name, " high (\u2265 ", threshold_label, ")")),
      legend.title = paste0(signature_name, " status"),
      palette = c("red", "blue"),
      title = paste("Kaplan-Meier Curve:", signature_name, "High vs Low", 
                    paste0("(", dataset, " data, ", threshold_label, " binary)")),
      tables.y.text = FALSE, 
      log.rank.weights = "1",
      ggtheme = theme_bw(),
      font.legend = font_legend,
      pval.size = risk_table_number_size,
      risk.table.fontsize = risk_table_number_size,
      align = TRUE, 
      tables.theme = theme_cleantable() +
        theme(
          axis.title.x = element_blank(), 
          axis.title.y = element_blank(), 
          plot.title = element_blank(),
          axis.text.x = element_text(size = font, face = "plain", color = "black"),
          axis.text.y = element_text(size = font, face = "plain", color = "black")
        )
    )
    
    # Ensure main plot has black text and bold title
    plot$plot <- plot$plot + 
      theme(
        plot.title = element_text(size = font, face = "bold", color = "black"),
        axis.title.x = element_text(size = font, face = "plain", color = "black"),
        axis.title.y = element_text(size = font, face = "plain", color = "black"),
        axis.text.x = element_text(size = font, face = "plain", color = "black"),
        axis.text.y = element_text(size = font, face = "plain", color = "black")
      ) +
      annotate("text", x = 10, y = 0.125, label = p_label, parse = TRUE, size = font / 2.833)
    
    
    return(plot)
  }
  
  # Create binary groups with 0.25 cutoff
  binary_groups <- create_binary_groups(0.25, "0.25")
  binary_df <- binary_groups$combined_df
  
  # Perform survival analysis for 0.25 cutoff
  binary_hazard_ratio <- coxph(Surv(overall_survival_months, status_numeric) ~ group, data = binary_df)
  binary_surv_fit <- survfit(Surv(overall_survival_months, status_numeric) ~ group, data = binary_df)
  binary_surv_diff <- survdiff(Surv(overall_survival_months, status_numeric) ~ group, data = binary_df)
  binary_plot <- create_survival_plot(binary_surv_fit, binary_surv_diff, binary_df, "0.25", dataset)
  
  # Determine optimal cutoff
  optimal_cutoff <- switch(signature_name,
                           "SBS4" = SBS4_cutoff,
                           "SBS13" = SBS13_cutoff,
                           cutoff)
  
  # Initialize results with base components
  results <- list(
    plot_0.25 = binary_plot,
    hazard_0.25 = binary_hazard_ratio,
    signature_Low = binary_groups$low_group,
    signature_High = binary_groups$high_group
  )
  
  # Add optimal cutoff analysis if conditions are met
  if (dataset %in% c("paper", "external", "all") && !is.null(optimal_cutoff)) {
    # Create optimal binary groups
    optimal_groups <- create_binary_groups(optimal_cutoff, optimal_cutoff)
    optimal_df <- optimal_groups$combined_df
    
    # Perform survival analysis for optimal cutoff
    optimal_hazard_ratio <- coxph(Surv(overall_survival_months, status_numeric) ~ group, data = optimal_df)
    optimal_surv_fit <- survfit(Surv(overall_survival_months, status_numeric) ~ group, data = optimal_df)
    optimal_surv_diff <- survdiff(Surv(overall_survival_months, status_numeric) ~ group, data = optimal_df)
    optimal_plot <- create_survival_plot(optimal_surv_fit, optimal_surv_diff, optimal_df, optimal_cutoff, dataset)
    
    # Add optimal results to the results list
    results <- c(results, list(
      plot_optimal = optimal_plot,
      hazard_optimal = optimal_hazard_ratio,
      optimal_signature_Low = optimal_groups$low_group,
      optimal_signature_High = optimal_groups$high_group
    ))
  }
  
  return(results)
}

analyze_signature_continuous <- function(signature_name, biological_weight_table, survival_data) {
  
  font <- 12
  font_legend <- 9
  risk_table_number_size <- font / 2.5
  
  dt <- merge(survival_data, biological_weight_table, by = "Unique_Patient_Identifier")
  
  if (!signature_name %in% colnames(dt)) {
    stop(paste("Column", signature_name, "not found in merged data"))
  }
  
  formula_str <- reformulate(signature_name, "Surv(overall_survival_months, status_numeric)")
  model <- coxph(formula_str, data = dt, x = TRUE)
  
  create_survival_plot_continuous <- plot_surv_area(
    time = "overall_survival_months",
    status = "status_numeric",
    variable = signature_name,
    data = dt,
    model = model,
    discrete = FALSE,
    start_color = "red",
    end_color = "blue",
    xlab = "Time (months)",
    ylab = "Overall survival",
    label_digits = 2,
    legend.title = element_blank(),  # uncomment for no legend title
    # legend.title = paste0(signature_name, " signature activity")
  ) +
    theme_bw() +
    theme(
      plot.title = element_blank(),
      axis.title.x = element_text(size = font, face = "plain", color = "black"),
      axis.title.y = element_text(size = font, face = "plain", color = "black"),
      axis.text.x = element_text(size = font, face = "plain", color = "black"),
      axis.text.y = element_text(size = font, face = "plain", color = "black"),
      legend.text = element_text(size = font_legend, color = "black"),
      legend.position = "top",
      legend.key.height = unit(0.3, "cm")
    ) +
    scale_x_continuous(limits = c(0, 80), breaks = seq(0, 80, 20)) +
    scale_y_continuous(limits = c(0, 1), breaks = seq(0, 1, 0.25))
  
  
  return(list(model_continuous = model, plot_continuous = create_survival_plot_continuous))
}

analyze_signature_complete <- function(signature_name, datasets = c("external", "paper", "all"), cutoff = NULL) {
  
  results <- list()
  
  dataset_config <- list(
    external = list(
      biological_weight_table = biological_weights_external,
      survival_data = survival_data_external
    ),
    paper = list(
      biological_weight_table = biological_weights_paper,
      survival_data = survival_data_paper
    ),
    all = list(
      biological_weight_table = biological_weights_all,
      survival_data = survival_data_all
    )
  )
  
  for (dataset in datasets) {
    cat(paste("Analyzing", signature_name, "for", dataset, "dataset...\n"))
    
    sig_analysis <- dataset_config[[dataset]]$biological_weight_table
    surv_data <- dataset_config[[dataset]]$survival_data
    
    tryCatch({
      sig_result_binary <- analyze_signature(
        signature_name = signature_name, 
        dataset = dataset, 
        biological_weight_table = sig_analysis, 
        survival_data = surv_data, 
        cutoff = cutoff
      )
      
      sig_result_continuous <- analyze_signature_continuous(
        signature_name = signature_name, 
        biological_weight_table = sig_analysis,
        survival_data = surv_data
      )
      
      results[[dataset]] <- list(
        binary = sig_result_binary,
        continuous = sig_result_continuous
      )
      
      cat(paste("Successfully analyzed", signature_name, "for", dataset, "dataset\n"))
      
    }, error = function(e) {
      cat(paste("Skipping", dataset, "due to error:", e$message, "\n"))
      results[[dataset]] <- NULL
    })
  }
  
  return(results)
}

create_diverging_signature_plot <- function(Li_data, ref_based_data) {
  
  font <- 11
  scale <- 7.25
  
  # LEFT SIDE: Signature groupings for shared signatures
  left_signature_groupings <- list(
    "Unknown, clock-like (5)" = "SBS5",
    "APOBEC (13)" = c("SBS13"),
    "Tobacco (4)" = c("SBS4")
  )
  left_other_label <- "Non-actionable and unknown"
  
  # LEFT SIDE: Cannataro colors
  left_cannataro_colors <- c(
    "Unknown, clock-like (5)" = "gray60",
    "APOBEC (13)" = "#7570b3",
    "Tobacco (4)" = "#a6761d",
    "Non-actionable and unknown" = "black"
  )
  
  # RIGHT SIDE: Signature groupings for ref_based only signatures
  right_signature_groupings <- list(
    "Deamination, clock-like (1)" = "SBS1",
    "APOBEC (2)" = c("SBS2"),
    "Defective homologous recombination (3)" = "SBS3",
    "Tobacco (29)" = c("SBS29"),
    "UV light (7a–d, 38)" = c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"),
    "Treatment (31, 32, 35, 86, 87)" = c("SBS11", "SBS31", "SBS32", "SBS35", "SBS86", "SBS87", "SBS90"),
    "Mutagenic chemical exposure (24)" = c("SBS22", "SBS24", "SBS42", "SBS88"),
    "Alcohol-associated (16)" = "SBS16"
  )
  right_other_label <- "Non-actionable and unknown"
  
  # RIGHT SIDE: Cannataro colors
  right_cannataro_colors <- c(
    "Deamination, clock-like (1)" = "gray40",
    "APOBEC (2)" = "#7570b3",
    "Defective homologous recombination (3)" = "#e7298a",
    "Tobacco (29)" = "#a6761d",
    "UV light (7a–d, 38)" = "#e6ab02",
    "Treatment (31, 32, 35, 86, 87)" = "#1b9e77",
    "Mutagenic chemical exposure (24)" = "#66a61e",
    "Alcohol-associated (16)" = "#d95f02",
    "Non-actionable and unknown" = "black"
  )
  
  # Process Li_data (genome mutagenesis share)
  Li_data_dt <- as.data.table(Li_data)
  mean_Li_data <- Li_data_dt[, sapply(.SD, mean), .SDcols = names(Li_data_dt)[-c(1)]]
  df_Li <- data.table(type = 'Li_results', prop = mean_Li_data, name = names(mean_Li_data))
  
  # Process ref_based_data (genome mutagenesis share)
  ref_based_dt <- as.data.table(ref_based_data)
  not_zero <- ref_based_dt[, (sapply(.SD, function(x) any(x != 0))), .SDcols = names(ref_based_dt)[-c(1:4)]]
  not_zero <- names(which(not_zero == TRUE))
  signature_weights <- ref_based_dt[, c("Unique_Patient_Identifier", ..not_zero)]
  signature_names <- names(signature_weights)[-c(1)]
  mean_signature_weights <- signature_weights[, sapply(.SD, mean), .SDcols = signature_names]
  df_ref <- data.table(type = 'ref_based', prop = mean_signature_weights, name = names(mean_signature_weights))
  
  # Combine data
  df_final <- rbind(df_Li, df_ref)
  
  # Determine which signatures are in both datasets vs only in ref_based
  li_signatures <- unique(df_Li$name)
  ref_signatures <- unique(df_ref$name)
  shared_signatures <- intersect(li_signatures, ref_signatures)
  ref_only_signatures <- setdiff(ref_signatures, li_signatures)
  
  # LEFT PANEL: Shared signatures (both Li_results and ref_based)
  left_dt <- df_final[name %in% shared_signatures]
  
  # Add signature groupings to other signatures for LEFT panel
  left_other_signatures <- setdiff(left_dt$name, unlist(left_signature_groupings))
  left_signature_groupings[[left_other_label]] <- left_other_signatures
  
  # Unwind list to get a table matching signatures to their labels for LEFT panel
  left_signature_labels <- rbindlist(lapply(1:length(left_signature_groupings), 
                                            function(x) data.table(name = left_signature_groupings[[x]], 
                                                                   label = names(left_signature_groupings)[x])))
  
  # Create group column and assign signature_groupings for LEFT panel
  left_dt[left_signature_labels, group := label, on = 'name']
  # Assign any unmatched signatures to the "other" category
  left_dt[is.na(group), group := left_other_label]
  
  # Order signature groups by mean effect share for LEFT panel
  left_summed_shares <- left_dt[type == 'ref_based' & group != left_other_label, .(summed_share = sum(prop)), by = 'group']
  left_mean_share_order <- left_summed_shares[order(summed_share), group]
  left_final_order <- c(left_other_label, left_mean_share_order)
  left_dt$group <- factor(left_dt$group, levels = left_final_order)
  
  left_dt[, type := factor(type, levels = c("Li_results", "ref_based"))]
  
  p_left <- ggplot(data = left_dt) +
    geom_col(
      mapping = aes(
        x = prop,
        y = type,
        fill = group#,
        #weight = prop
      ),
      position = position_stack(reverse=FALSE),
      width = .8,
      color = 'black'
    ) + labs(x="Proportion of TMB", y="") +
    scale_x_reverse(
      limits = c(1, 0),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) format(abs(x), nsmall = 2),
      expand = expansion(mult = c(0.05, 0))  # 5% expansion on left (1.00), no expansion on right (0.00)
    ) +
    theme_bw(base_size = font) +
    theme(
      axis.title = element_blank(),
      axis.title.x = element_text(color = "black", size = 12, face="plain"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 2, 0)
    )
  
  # RIGHT PANEL: ref_based only signatures
  right_dt <- df_final[name %in% ref_only_signatures]
  
  # Add signature groupings to other signatures for RIGHT panel
  right_other_signatures <- setdiff(right_dt$name, unlist(right_signature_groupings))
  right_signature_groupings[[right_other_label]] <- right_other_signatures
  
  # Unwind list to get a table matching signatures to their labels for RIGHT panel
  right_signature_labels <- rbindlist(lapply(1:length(right_signature_groupings), 
                                             function(x) data.table(name = right_signature_groupings[[x]], 
                                                                    label = names(right_signature_groupings)[x])))
  
  # Add zero entries for Li_results to maintain structure
  if(nrow(right_dt) > 0) {
    zero_entries <- data.table(
      type = 'Li_results',
      prop = 0,
      name = unique(right_dt$name)
    )
    right_dt_extended <- rbind(right_dt, zero_entries)
    right_dt_extended[, type := factor(type, levels = c("Li_results", "ref_based"))]
    
    # Create group column and assign signature_groupings for RIGHT panel
    right_dt_extended[right_signature_labels, group := label, on = 'name']
    # Assign any unmatched signatures to the "other" category
    right_dt_extended[is.na(group), group := right_other_label]
    
    # Order signature groups by mean effect share for RIGHT panel
    right_summed_shares <- right_dt_extended[type == 'ref_based' & group != right_other_label, .(summed_share = sum(prop)), by = 'group']
    right_mean_share_order <- right_summed_shares[order(summed_share), group]
    right_final_order <- c(right_other_label, right_mean_share_order)
    right_dt_extended$group <- factor(right_dt_extended$group, levels = right_final_order)
  } else {
    right_dt_extended <- data.table(type = character(), prop = numeric(), name = character(), group = character())
  }
  
  p_right <- ggplot(data = right_dt_extended) +
    geom_col(
      mapping = aes(
        x = prop,
        y = type,
        fill = group),
      position = position_stack(reverse = FALSE),
      width = .8,
      color = 'black'
    ) + labs(x="Proportion of TMB", y="") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) format(abs(x), nsmall = 2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_bw(base_size = font) +
    theme(
      axis.title = element_blank(),
      axis.title.x = element_text(color = "black", size = 12, face="plain"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 2, -3)
    )
  
  # LABEL PANEL
  label_dt <- data.table(type = c(0, 0.7), label = c("atop(italic('De novo'), plain('extraction'))", "atop(plain('Reference'), plain('based'))"))
  
  p_label <- ggplot(label_dt) +
    geom_text(aes(x = 0, y = type, label = label),
              size = 12 / 2.8, 
              fontface = "plain", 
              hjust = 0.5, parse = TRUE) +
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    scale_y_continuous(limits = c(-1, 1)) +  # Add y scale to allow positioning
    theme_void() +
    theme(
      # plot.margin = margin(5, 2, 5, 2),
      plot.margin = margin(0, 0, 0, 2),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  
  # Set legend name
  sig_legend_name <- ''
  
  # Order groups for left panel legend
  if(nrow(left_dt) > 0) {
    left_group_totals <- left_dt[, .(total_proportion = sum(prop)), by = group]
    left_group_order <- left_group_totals[order(-total_proportion), group]
  } else {
    left_group_order <- unique(left_dt$group)
  }
  
  # Order groups for right panel legend
  if(nrow(right_dt_extended) > 0) {
    right_group_totals <- right_dt_extended[, .(total_proportion = sum(prop)), by = group]
    right_group_order <- right_group_totals[order(-total_proportion), group]
  } else {
    right_group_order <- unique(right_dt_extended$group)
  }
  
  # Apply legend styling to left plot
  p_left <- p_left +
    scale_fill_manual(name = sig_legend_name, values = left_cannataro_colors,
                      breaks = left_final_order) +
    guides(fill = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(color = "black", size = scale, face="plain"), legend.key.size = unit(scale/10, "lines"))
  
  # Apply legend styling to right plot
  p_right <- p_right +
    scale_fill_manual(name = sig_legend_name, values = right_cannataro_colors,
                      breaks = right_final_order) +
    guides(fill = guide_legend(reverse = TRUE, nrow = 4, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(color = "black", size = scale, face="plain"), legend.key.size = unit(scale/10, "lines"))
  # https://github.com/wilkelab/cowplot/issues/202#issuecomment-2550098822
  get_legend2 <- function(plot, legend = NULL) {
    gt <- ggplotGrob(plot)
    pattern <- "guide-box"
    if (!is.null(legend)) {
      pattern <- paste0(pattern, "-", legend)
    }
    indices <- grep(pattern, gt$layout$name)
    not_empty <- !vapply(
      gt$grobs[indices], 
      inherits, what = "zeroGrob", 
      FUN.VALUE = logical(1)
    )
    indices <- indices[not_empty]
    if (length(indices) > 0) {
      return(gt$grobs[[indices[1]]])
    }
    return(NULL)
  }
  legend_left <- get_legend2(p_left, legend = "bottom")
  legend_right <- get_legend2(p_right, legend = "bottom")
  panel_1 <- plot_grid(NULL, legend_left, legend_right, nrow = 1, ncol = 3, rel_widths = c(0.75, 3.5, 4.8))
  p_right <- p_right + theme(legend.position = "none")
  p_left <- p_left + theme(legend.position = "none")
  # #combine with plot_grid
  panel_2<- plot_grid(p_label, p_left, p_right, NULL, ncol =4, nrow = 1, rel_widths = c(1.2,4.2,4.2,1.2))
  combined_plots <- plot_grid(panel_1, panel_2, nrow=2, ncol=1, rel_heights = c(7,12))
  
  return(combined_plots)
}


create_diverging_signature_plot_with_cancer_effect <- function(Li_data, ref_based_data, effect_data) {
  
  font <- 11
  scale <- 7.25
  
  # LEFT SIDE: Signature groupings for shared signatures
  left_signature_groupings <- list(
    "Unknown, clock-like (5)" = "SBS5",
    "APOBEC (13)" = c("SBS13"),
    "Tobacco (4)" = c("SBS4")
  )
  left_other_label <- "Non-actionable and unknown"
  
  # LEFT SIDE: Cannataro colors
  left_cannataro_colors <- c(
    "Unknown, clock-like (5)" = "gray60",
    "APOBEC (13)" = "#7570b3",
    "Tobacco (4)" = "#a6761d",
    "Non-actionable and unknown" = "black"
  )
  
  # RIGHT SIDE: Signature groupings for ref_based only signatures
  right_signature_groupings <- list(
    "Deamination, clock-like (1)" = "SBS1",
    "APOBEC (2)" = c("SBS2"),
    "Defective homologous recombination (3)" = "SBS3",
    "Tobacco (29)" = c("SBS29"),
    "UV light (7a–d, 38)" = c("SBS7a", "SBS7b", "SBS7c", "SBS7d", "SBS38"),
    "Treatment (31, 32, 35, 86, 87)" = c("SBS11", "SBS31", "SBS32", "SBS35", "SBS86", "SBS87", "SBS90"),
    "Mutagenic chemical exposure (24)" = c("SBS22", "SBS24", "SBS42", "SBS88"),
    "Alcohol-associated (16)" = "SBS16"
  )
  right_other_label <- "Non-actionable and unknown"
  
  # RIGHT SIDE: Cannataro colors
  right_cannataro_colors <- c(
    "Deamination, clock-like (1)" = "gray40",
    "APOBEC (2)" = "#7570b3",
    "Defective homologous recombination (3)" = "#e7298a",
    "Tobacco (29)" = "#a6761d",
    "UV light (7a–d, 38)" = "#e6ab02",
    "Treatment (31, 32, 35, 86, 87)" = "#1b9e77",
    "Mutagenic chemical exposure (24)" = "#66a61e",
    "Alcohol-associated (16)" = "#d95f02",
    "Non-actionable and unknown" = "black"
  )
  
  # Process Li_data (genome mutagenesis share)
  Li_data_dt <- as.data.table(Li_data)
  mean_Li_data <- Li_data_dt[, sapply(.SD, mean), .SDcols = names(Li_data_dt)[-c(1)]]
  df_Li <- data.table(type = 'Li_results', prop = mean_Li_data, name = names(mean_Li_data))
  
  # Process ref_based_data (genome mutagenesis share)
  ref_based_dt <- as.data.table(ref_based_data)
  not_zero <- ref_based_dt[, (sapply(.SD, function(x) any(x != 0))), .SDcols = names(ref_based_dt)[-c(1:4)]]
  not_zero <- names(which(not_zero == TRUE))
  signature_weights <- ref_based_dt[, c("Unique_Patient_Identifier", ..not_zero)]
  signature_names <- names(signature_weights)[-c(1)]
  mean_signature_weights <- signature_weights[, sapply(.SD, mean), .SDcols = signature_names]
  df_ref <- data.table(type = 'ref_based', prop = mean_signature_weights, name = names(mean_signature_weights))
  
  # Combine data
  df_final <- rbind(df_ref, df_Li)
  
  # Determine which signatures are in both datasets vs only in ref_based
  li_signatures <- unique(df_Li$name)
  ref_signatures <- unique(df_ref$name)
  shared_signatures <- intersect(li_signatures, ref_signatures)
  ref_only_signatures <- setdiff(ref_signatures, li_signatures)
  
  # LEFT PANEL: Shared signatures (both Li_results and ref_based)
  left_dt <- df_final[name %in% shared_signatures]
  
  # Add signature groupings to other signatures for LEFT panel
  left_other_signatures <- setdiff(left_dt$name, unlist(left_signature_groupings))
  left_signature_groupings[[left_other_label]] <- left_other_signatures
  
  # Unwind list to get a table matching signatures to their labels for LEFT panel
  left_signature_labels <- rbindlist(lapply(1:length(left_signature_groupings), 
                                            function(x) data.table(name = left_signature_groupings[[x]], 
                                                                   label = names(left_signature_groupings)[x])))
  
  # Create group column and assign signature_groupings for LEFT panel
  left_dt[left_signature_labels, group := label, on = 'name']
  # Assign any unmatched signatures to the "other" category
  left_dt[is.na(group), group := left_other_label]
  
  # Order signature groups by mean effect share for LEFT panel
  left_summed_shares <- left_dt[type == 'ref_based' & group != left_other_label, .(summed_share = sum(prop)), by = 'group']
  left_mean_share_order <- left_summed_shares[order(summed_share), group]
  left_final_order <- c(left_other_label, left_mean_share_order)
  left_dt$group <- factor(left_dt$group, levels = left_final_order)
  
  left_dt[, type := factor(type, levels = c("ref_based", "Li_results"))]
  
  p_left <- ggplot(data = left_dt) +
    geom_col(
      mapping = aes(
        x = prop,
        y = type,
        fill = group#,
        #weight = prop
      ),
      position = position_stack(reverse=FALSE),
      width = .8,
      color = 'black'
    ) + labs(x="Proportion of TMB", y="") +
    scale_x_reverse(
      limits = c(1, 0),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) format(abs(x), nsmall = 2),
      expand = expansion(mult = c(0.05, 0))  # 5% expansion on left (1.00), no expansion on right (0.00)
    ) +
    theme_bw(base_size = font) +
    theme(
      axis.title = element_blank(),
      axis.title.x = element_text(color = "black", size = 12, face="plain"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 2, 0)
    )
  
  # RIGHT PANEL: ref_based only signatures
  right_dt <- df_final[name %in% ref_only_signatures]
  
  # Add signature groupings to other signatures for RIGHT panel
  right_other_signatures <- setdiff(right_dt$name, unlist(right_signature_groupings))
  right_signature_groupings[[right_other_label]] <- right_other_signatures
  
  # Unwind list to get a table matching signatures to their labels for RIGHT panel
  right_signature_labels <- rbindlist(lapply(1:length(right_signature_groupings), 
                                             function(x) data.table(name = right_signature_groupings[[x]], 
                                                                    label = names(right_signature_groupings)[x])))
  
  # Add zero entries for Li_results to maintain structure
  if(nrow(right_dt) > 0) {
    zero_entries <- data.table(
      type = 'Li_results',
      prop = 0,
      name = unique(right_dt$name)
    )
    right_dt_extended <- rbind(right_dt, zero_entries)
    right_dt_extended[, type := factor(type, levels = c("ref_based", "Li_results"))]
    
    # Create group column and assign signature_groupings for RIGHT panel
    right_dt_extended[right_signature_labels, group := label, on = 'name']
    # Assign any unmatched signatures to the "other" category
    right_dt_extended[is.na(group), group := right_other_label]
    
    # Order signature groups by mean effect share for RIGHT panel
    right_summed_shares <- right_dt_extended[type == 'ref_based' & group != right_other_label, .(summed_share = sum(prop)), by = 'group']
    right_mean_share_order <- right_summed_shares[order(summed_share), group]
    right_final_order <- c(right_other_label, right_mean_share_order)
    right_dt_extended$group <- factor(right_dt_extended$group, levels = right_final_order)
  } else {
    right_dt_extended <- data.table(type = character(), prop = numeric(), name = character(), group = character())
  }
  
  p_right <- ggplot(data = right_dt_extended) +
    geom_col(
      mapping = aes(
        x = prop,
        y = type,
        fill = group),
      position = position_stack(reverse = FALSE),
      width = .8,
      color = 'black'
    ) + labs(x="Proportion of TMB", y="") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) format(abs(x), nsmall = 2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_bw(base_size = font) +
    theme(
      axis.title = element_blank(),
      axis.title.x = element_text(color = "black", size = 12, face="plain"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 2, -3)
    )
  
  # LABEL PANEL
  label_dt <- data.table(type = c(0, 0.7), label = c("atop(plain('Reference'), plain('based'))", "atop(italic('De novo'), plain('extraction'))"))
  
  p_label <- ggplot(label_dt) +
    geom_text(aes(x = 0, y = type, label = label),
              size = 12 / 2.8, 
              fontface = "plain", 
              hjust = 0.5, parse = TRUE) +
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    scale_y_continuous(limits = c(-1, 1)) +  # Add y scale to allow positioning
    theme_void() +
    theme(
      # plot.margin = margin(5, 2, 5, 2),
      plot.margin = margin(0, 0, 0, 2),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  
  # Set legend name
  sig_legend_name <- ''
  
  # Order groups for left panel legend
  if(nrow(left_dt) > 0) {
    left_group_totals <- left_dt[, .(total_proportion = sum(prop)), by = group]
    left_group_order <- left_group_totals[order(-total_proportion), group]
  } else {
    left_group_order <- unique(left_dt$group)
  }
  
  # Order groups for right panel legend
  if(nrow(right_dt_extended) > 0) {
    right_group_totals <- right_dt_extended[, .(total_proportion = sum(prop)), by = group]
    right_group_order <- right_group_totals[order(-total_proportion), group]
  } else {
    right_group_order <- unique(right_dt_extended$group)
  }
  
  # Apply legend styling to left plot
  p_left <- p_left +
    scale_fill_manual(name = sig_legend_name, values = left_cannataro_colors,
                      breaks = left_final_order) +
    guides(fill = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(color = "black", size = scale, face="plain"), legend.key.size = unit(scale/10, "lines"))
  
  # Apply legend styling to right plot
  p_right <- p_right +
    scale_fill_manual(name = sig_legend_name, values = right_cannataro_colors,
                      breaks = right_final_order) +
    guides(fill = guide_legend(reverse = TRUE, nrow = 4, byrow = TRUE)) +
    theme(legend.position = "bottom", legend.text = element_text(color = "black", size = scale, face="plain"), legend.key.size = unit(scale/10, "lines"))
  # https://github.com/wilkelab/cowplot/issues/202#issuecomment-2550098822
  get_legend2 <- function(plot, legend = NULL) {
    gt <- ggplotGrob(plot)
    pattern <- "guide-box"
    if (!is.null(legend)) {
      pattern <- paste0(pattern, "-", legend)
    }
    indices <- grep(pattern, gt$layout$name)
    not_empty <- !vapply(
      gt$grobs[indices], 
      inherits, what = "zeroGrob", 
      FUN.VALUE = logical(1)
    )
    indices <- indices[not_empty]
    if (length(indices) > 0) {
      return(gt$grobs[[indices[1]]])
    }
    return(NULL)
  }
  legend_left <- get_legend2(p_left, legend = "bottom")
  legend_right <- get_legend2(p_right, legend = "bottom")
  panel_1 <- plot_grid(NULL, legend_left, legend_right, nrow = 1, ncol = 3, rel_widths = c(0.75, 3.5, 4.8))
  p_right <- p_right + theme(legend.position = "none")
  p_left <- p_left + theme(legend.position = "none")
  # #combine with plot_grid
  panel_2<- plot_grid(p_label, p_left, p_right, NULL, ncol =4, nrow = 1, rel_widths = c(1.2,4.2,4.2,1.2))
  
  # Process effect_data (effect share)
  # Check if effect_data is a named vector (average shares) or data.table/data.frame
  if(is.vector(effect_data) && !is.null(names(effect_data))) {
    # effect_data is already the average effect shares vector
    mean_effect_weights <- effect_data
    df_effect <- data.table(type = 'effect', prop = mean_effect_weights, name = names(mean_effect_weights))
  } else {
    # effect_data is a data.table/data.frame, process as before
    effect_dt <- as.data.table(effect_data)
    not_zero <- effect_dt[, (sapply(.SD, function(x) any(x != 0))), .SDcols = names(effect_dt)[-c(1:4)]]
    not_zero <- names(which(not_zero == TRUE))
    signature_weights <- effect_dt[, c("Unique_Patient_Identifier", ..not_zero)]
    signature_names <- names(signature_weights)[-c(1)]
    mean_effect_weights <- signature_weights[, sapply(.SD, mean), .SDcols = signature_names]
    df_effect <- data.table(type = 'effect', prop = mean_effect_weights, name = names(mean_effect_weights))
  }
  
  # LEFT PANEL: 5 identified signatures
  left_dt <- df_effect[name %in% shared_signatures]
  
  # Add signature groupings to other signatures for LEFT panel
  left_other_signatures <- setdiff(left_dt$name, unlist(left_signature_groupings))
  left_signature_groupings[[left_other_label]] <- left_other_signatures
  
  # Create group column and assign signature_groupings for LEFT panel
  left_dt[left_signature_labels, group := label, on = 'name']
  # Assign any unmatched signatures to the "other" category
  left_dt[is.na(group), group := left_other_label]
  
  # Order signature groups by mean effect share for LEFT panel
  left_summed_shares <- left_dt[type == 'effect' & group != left_other_label, .(summed_share = sum(prop)), by = 'group']
  left_mean_share_order <- left_summed_shares[order(summed_share), group]
  left_final_order <- c(left_other_label, left_mean_share_order)
  left_dt$group <- factor(left_dt$group, levels = left_final_order)
  
  left_dt[, type := factor(type, levels = "effect")]
  
  p_left_effect <- ggplot(data = left_dt) +
    geom_col(
      mapping = aes(
        x = prop,
        y = type,
        fill = group
      ),
      position = position_stack(reverse=FALSE),
      width = .8,
      color = 'black'
    ) + labs(x="Proportion of oncogenic effect", y="") +
    scale_x_reverse(
      limits = c(1, 0),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) format(abs(x), nsmall = 2),
      expand = expansion(mult = c(0.05, 0))  # 5% expansion on left (1.00), no expansion on right (0.00)
    ) +
    theme_bw(base_size = font) +
    theme(
      axis.title = element_blank(),
      axis.title.x = element_text(color = "black", size = 12, face="plain"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 2, 0)
    )
  
  # RIGHT PANEL: ref_based only signatures for EFFECT data
  right_dt <- df_effect[name %in% ref_only_signatures]
  
  # Add signature groupings to other signatures for RIGHT panel
  right_other_signatures <- setdiff(right_dt$name, unlist(right_signature_groupings))
  right_signature_groupings[[right_other_label]] <- right_other_signatures
  
  if(nrow(right_dt) > 0) {
    # Create group column and assign signature_groupings for RIGHT panel
    right_dt[right_signature_labels, group := label, on = 'name']
    # Assign any unmatched signatures to the "other" category
    right_dt[is.na(group), group := right_other_label]
    
    # Order signature groups by mean effect share for RIGHT panel
    right_summed_shares <- right_dt[group != right_other_label, .(summed_share = sum(prop)), by = 'group']
    right_mean_share_order <- right_summed_shares[order(summed_share), group]
    right_final_order <- c(right_other_label, right_mean_share_order)
    right_dt$group <- factor(right_dt$group, levels = right_final_order)
    
    right_dt[, type := factor(type, levels = "effect")]
  } else {
    right_dt <- data.table(type = character(), prop = numeric(), name = character(), group = character())
  }
  
  p_right_effect <- ggplot(data = right_dt) +
    geom_col(
      mapping = aes(
        x = prop,
        y = type,
        fill = group),
      position = position_stack(reverse = FALSE),
      width = .8,
      color = 'black'
    ) + labs(x="Proportion of oncogenic effect", y="") +
    scale_x_continuous(
      limits = c(0, 1),
      breaks = seq(0, 1, by = 0.25),
      labels = function(x) format(abs(x), nsmall = 2),
      expand = expansion(mult = c(0, 0.05))
    ) +
    theme_bw(base_size = font) +
    theme(
      axis.title = element_blank(),
      axis.title.x = element_text(color = "black", size = 12, face="plain"),
      axis.text.y = element_blank(),
      axis.ticks.y = element_blank(),
      axis.text.x = element_text(size = 12),
      panel.grid = element_blank(),
      legend.position = "none",
      plot.margin = margin(0, 0, 2, -3)
    )
  
  # LABEL PANEL for effect
  label_dt_effect <- data.table(type = 0.4, label = "atop(plain('Reference'), plain('based'))")
  
  p_label_effect <- ggplot(label_dt_effect) +
    geom_text(aes(x = 0, y = type, label = label),
              size = 12 / 2.8, 
              fontface = "plain", 
              hjust = 0.5, parse = TRUE) +  # Added parse = TRUE for atop()
    scale_x_continuous(limits = c(-0.5, 0.5)) +
    scale_y_continuous(limits = c(-1, 1)) +
    theme_void() +
    theme(
      plot.margin = margin(0, 0, 0, 2),
      panel.background = element_blank(),
      plot.background = element_blank()
    )
  
  # Set legend name
  sig_legend_name <- ''
  
  # Order groups for left panel legend
  if(nrow(left_dt) > 0) {
    left_group_totals <- left_dt[, .(total_proportion = sum(prop)), by = group]
    left_group_order <- left_group_totals[order(-total_proportion), group]
  } else {
    left_group_order <- unique(left_dt$group)
  }
  
  # Order groups for right panel legend
  if(nrow(right_dt_extended) > 0) {
    right_group_totals <- right_dt_extended[, .(total_proportion = sum(prop)), by = group]
    right_group_order <- right_group_totals[order(-total_proportion), group]
  } else {
    right_group_order <- unique(right_dt_extended$group)
  }
  
  # Apply legend styling to left plot
  p_left_effect <- p_left_effect +
    scale_fill_manual(name = sig_legend_name, values = left_cannataro_colors,
                      breaks = left_final_order) +
    guides(fill = guide_legend(reverse = TRUE, nrow = 2, byrow = TRUE)) +
    theme(legend.position = "none", legend.text = element_text(color = "black", size = scale, face="plain"), legend.key.size = unit(scale/10, "lines"))
  
  # Apply legend styling to right plot
  p_right_effect <- p_right_effect +
    scale_fill_manual(name = sig_legend_name, values = right_cannataro_colors,
                      breaks = right_final_order) +
    guides(fill = guide_legend(reverse = TRUE, nrow = 4, byrow = TRUE)) +
    theme(legend.position = "none", legend.text = element_text(color = "black", size = scale, face="plain"), legend.key.size = unit(scale/10, "lines"))
  
  panel_3<- plot_grid(p_label_effect, p_left_effect, p_right_effect, NULL, ncol =4, nrow = 1, rel_widths = c(1.2,4.2,4.2,1.2))
  
  combined_plots <- plot_grid(panel_1, panel_2, panel_3, nrow=3, ncol=1, rel_heights = c(7,12, 7))
  
  return(combined_plots)
}

# ----------------Gene Level Analysis ------------------------

# Define known gene classifications (Specified by COSMIC or OncoKB) for compound variant analysis

# Zhou_et_al characterized 9610 genes;
genes_info <- read.xlsx(Zhou_et_al_source_data, sheet = 7, startRow = 1, cols = 1:8) %>%
  filter(Role != "Passenger") %>%
  select(gene_name = gene, cancer_role = Role) %>%
  mutate(cancer_role = if_else(cancer_role == "Other", "TSG", cancer_role)) %>%
  distinct(gene_name, .keep_all = TRUE)

# Function to filter variants based on gene name and cancer role in order to run compound variant analysis
filter_variants_by_gene <- function(gene_name, cancer_role, variants_table = variants) {
  
  # Apply filtering based on cancer role
  if (cancer_role == "TSG") {
    # For tumor suppressor genes: include variants with maf_prevalence > 1 OR nonsense mutations, exclude intergenic
    filtered_variants <- variants_table[
      gene == gene_name & 
        (maf_prevalence > 1 | 
           (aa_ref != "STOP" & aa_alt == "STOP") | 
           (aa_ref == "STOP" & aa_alt != "STOP")) & 
        intergenic == FALSE
    ]
  } else {
    # For oncogenes: only include variants with maf_prevalence > 1
    filtered_variants <- variants_table[
      gene == gene_name & 
        maf_prevalence > 1
    ]
  }
  
  return(filtered_variants)
}

# Function to process all genes at once for compound variant analysis
process_all_genes <- function(top_genes_for_compound_variant, genes_info) {
  # Create internal_gene_info data.table
  internal_gene_info <- data.table(
    internal_gene_name = top_genes_for_compound_variant
  )
  
  # Join with genes_info on gene_name to get cancer_role
  # We'll do a left join and assign "TSG" if no match found
  internal_gene_info <- merge(
    internal_gene_info, 
    genes_info, 
    by.x = "internal_gene_name", 
    by.y = "gene_name", 
    all.x = TRUE
  )
  
  # Replace NA cancer_role with "TSG"
  internal_gene_info[is.na(cancer_role), cancer_role := "TSG"]
  
  # Initialize a vector to hold all top variants
  all_top_variants_list <- c()
  
  # Iterate over each gene/classification pair
  for (i in seq_len(nrow(internal_gene_info))) {
    gene <- as.character(internal_gene_info$internal_gene_name[i])
    role <- as.character(internal_gene_info$cancer_role[i])
    
    # Call the filter_variants_by_gene function
    top_variants <- filter_variants_by_gene(gene_name = gene, cancer_role = role, variants_table = variants)
    
    # Save the filtered variants data.table to the list
    all_top_variants_list[[i]] <- top_variants
  }
  
  # Combine all data.tables into one; handle empty list case safely
  if (length(all_top_variants_list) == 0) {
    return(data.table())  # empty data.table
  } else {
    combined_variants <- rbindlist(all_top_variants_list, use.names = TRUE, fill = TRUE)
    return(combined_variants)
  }
}

get_prognostic_significance <- function(input_df, top_x_genes, master_df, filter_by) {
  
  # Extract top genes based on selection intensity
  gene_list <- input_df %>%
    arrange(desc(selection_intensity)) %>%
    distinct(gene, .keep_all = TRUE) %>%
    slice_head(n = top_x_genes) %>%
    pull(gene)
  
  # Identify missing and matched genes
  missing_genes <- setdiff(gene_list, master_df$Gene)
  matched_genes <- intersect(gene_list, master_df$Gene)
  
  # Process data based on filter_by parameter (using 0.05 significance threshold)
  if (filter_by == "fdr") {
    # Genes significant in univariate (FDR < 0.05)
    sig_uni_genes <- master_df %>%
      filter(Univariate_Cox_fdr < 0.05) %>%
      select(Gene, Univariate_Cox_fdr)
    
    # Genes significant in multivariate (FDR < 0.05)
    sig_multi_genes <- master_df %>%
      filter(Multivariate_Cox_fdr < 0.05) %>%
      select(Gene, Multivariate_Cox_fdr)
    
    df <- master_df %>%
      filter(Gene %in% matched_genes) %>%
      select(Gene, Univariate_Cox_fdr, Multivariate_Cox_fdr) %>%
      mutate(
        univariate_significant_status = ifelse(Gene %in% sig_uni_genes$Gene, "yes", "no"),
        multivariate_significant_status = ifelse(Gene %in% sig_multi_genes$Gene, "yes", "no")
      ) %>%
      select(Gene, univariate_significant_status, multivariate_significant_status,
             univariate_fdr = Univariate_Cox_fdr,
             multivariate_fdr = Multivariate_Cox_fdr)
    
    df_uni_sig <- with(df, Gene[univariate_significant_status == "yes"])
    df_multi_sig <- with(df, Gene[multivariate_significant_status == "yes"])
  } else {
    # Genes significant in univariate (P < 0.05)
    sig_uni_genes <- master_df %>%
      filter(Univariate_Cox_P < 0.05) %>%
      select(Gene, Univariate_Cox_P)
    
    # Genes significant in multivariate (P < 0.05)
    sig_multi_genes <- master_df %>%
      filter(Multivariate_Cox_P < 0.05) %>%
      select(Gene, Multivariate_Cox_P)
    
    df <- master_df %>%
      filter(Gene %in% matched_genes) %>%
      select(Gene, Univariate_Cox_P, Multivariate_Cox_P) %>%
      mutate(
        univariate_significant_status = ifelse(Gene %in% sig_uni_genes$Gene, "yes", "no"),
        multivariate_significant_status = ifelse(Gene %in% sig_multi_genes$Gene, "yes", "no")
      ) %>%
      select(Gene, univariate_significant_status, multivariate_significant_status,
             univariate_P = Univariate_Cox_P,
             multivariate_P = Multivariate_Cox_P)
    
    df_uni_sig <- with(df, Gene[univariate_significant_status == "yes"])
    df_multi_sig <- with(df, Gene[multivariate_significant_status == "yes"])
  }
  
  # Create summary results
  df_results <- list(
    Genes_Uni_Sig = df_uni_sig,
    Genes_Multi_Sig = df_multi_sig,
    Genes_Missing = missing_genes,
    Proportion_Missing = length(missing_genes) / top_x_genes
  )
  
  # Return both detailed results and summary
  return(list(result = df, summary = df_results))
}

get_gene_prognosis_from_signature_attribution <- function(cesa, signature_attribution, signature) {
  df <- signature_attribution$mutational_sources$average_by_variant
  run_name <- paste(signature, "related variants")
  
  # Columns to compare against: all except first column and the signature itself
  cols_to_compare <- setdiff(names(df)[-1], signature)
  df$gene <- sub("_.*", "",df$variant_id)

  df <- df %>%
    rowwise() %>%
    filter(.data[[signature]] > max(c_across(all_of(cols_to_compare)), na.rm = TRUE)) %>%
    ungroup() %>%
    select("gene", "variant_id", all_of(signature))
  
  variants <- cesa$variants[cesa$variants$variant_id %in% df$variant_id,]
  cesa <- ces_variant(cesa = cesa, run_name = run_name, variants = variants)
  gene_prognosis_data <- get_prognostic_significance(cesa$selection[[run_name]], nrow(variants), cox_regression_df, filter_by = "P")
  return(list(
    filtered_df = df,
    gene_prognosis = gene_prognosis_data
  ))
}

