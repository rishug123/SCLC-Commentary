
##--------------- ANALYSIS FOR CESA_EXTERNAL (external data)-----------
# 
# Plot a linear regression of log(TMB) vs Mutations Attributed to SBS4 and SBS13 (external data)
external_TMB_linreg_SBS4 <- plot_TMB_linear_regression(biological_weights_external, "external", "SBS4")
external_TMB_linreg_SBS13 <- plot_TMB_linear_regression(biological_weights_external, "external", "SBS13")
external_TMB_linreg_SBS5 <- plot_TMB_linear_regression(biological_weights_external, "external", "SBS5")

# Plot Survival Curves 

# Signature Absence (<0.25 or <paper_bs_optimal_cutoff) vs Presence (>= 0.25 or >=paper_bs_optimal_cutoff)

SBS4_survival_0.25_external <- SBS4_results$external$binary$plot_0.25
SBS13_survival_0.25_external <- SBS13_results$external$binary$plot_0.25
SBS4_survival_optimal_external <- SBS4_results$external$binary$plot_optimal
SBS13_survival_optimal_external <- SBS13_results$external$binary$plot_optimal

##--------------- ANALYSIS FOR CESA_PAPER (paper data)-------
# Plot box plot of log(TMB) vs mutations attributed to SBS4 and SBS13
data_paper <- biological_weights_paper[, .(total_snvs, SBS3, SBS4, SBS5, SBS13, SBS24, SBS31, SBS87)]
data_paper <- data_paper[complete.cases(data_paper),]

SBS4_0.25 <- plot_TMB_binary_boxplot(data_paper, "SBS4", 0.25, "(Binary for Signature Presence: 0.25)")
SBS4_optimal <- plot_TMB_binary_boxplot(data_paper, "SBS4", SBS4_survival_optimal_cutoff_paper_mean, paste0("(Optimal Binary for Signature Presence: ",SBS4_survival_optimal_cutoff_paper_mean,")"))
SBS13_0.25 <- plot_TMB_binary_boxplot(data_paper, "SBS13", 0.25, "(Binary for Signature Presence: 0.25)", p_value = FALSE)
SBS13_optimal <- plot_TMB_binary_boxplot(data_paper, "SBS13", SBS13_survival_optimal_cutoff_paper_mean, paste0("(Optimal Binary for Signature Presence: ",SBS13_survival_optimal_cutoff_paper_mean,")"))

# Convert data to data.frames for ggplot2
sbs13_tmb_df <- data.frame(cutoff = SBS13_TMB_optimal_cutoff_paper[[1]])
sbs4_tmb_df  <- data.frame(cutoff = SBS4_TMB_optimal_cutoff_paper[[1]])
sbs13_surv_df <- data.frame(cutoff = SBS13_survival_optimal_cutoff_paper[[1]])
sbs4_surv_df  <- data.frame(cutoff = SBS4_survival_optimal_cutoff_paper[[1]])

font_size <- 12

# SBS4 x-axis range
sbs4_min <- 0
sbs4_max <- 0.6

# SBS13 x-axis range
sbs13_min <- 0
sbs13_max <- 0.15

# SBS4 TMB
sbs4_tmb_bs_distribution <- ggplot(sbs4_tmb_df, aes(x = cutoff)) +
  geom_histogram(fill = "#a6761d", color = "black", bins = 30) +
  geom_vline(aes(xintercept = SBS4_paper_optimal_TMB_cutoff_for_paper_non_bs$estimate, color = "Non-Bootstrapped Cutoff"), linewidth = 1) +
  geom_vline(aes(xintercept = SBS4_TMB_optimal_cutoff_paper_mean, color = "Bootstrap Mean"), size = 1) +
  scale_color_manual(name = NULL, values = c("Non-Bootstrapped Cutoff" = "red", "Bootstrap Mean" = "green")) +
  scale_x_continuous(limits = c(sbs4_min, sbs4_max)) +
  labs(title = "Bootstrap Distribution of Optimal SBS4 TMB Cutoff (Paper Data)", x = "SBS4 TMB cutoff value", y = "Frequency") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = font_size),
    axis.title = element_text(size = font_size)
  )

# SBS4 Survival
sbs4_surv_bs_distribution <- ggplot(sbs4_surv_df, aes(x = cutoff)) +
  geom_histogram(fill = "#a6761d", color = "black", bins = 30) +
  geom_vline(aes(xintercept = SBS4_paper_optimal_survival_cutoff_for_paper_non_bs$estimate, color = "Non-Bootstrapped Cutoff"), linewidth = 1) +
  geom_vline(aes(xintercept = SBS4_survival_optimal_cutoff_paper_mean, color = "Bootstrap Mean"), size = 1) +
  scale_color_manual(name = NULL, values = c("Non-Bootstrapped Cutoff" = "red", "Bootstrap Mean" = "green")) +
  scale_x_continuous(limits = c(sbs4_min, sbs4_max)) +
  labs(title = "Bootstrap Distribution of Optimal SBS4 Survival Cutoff (Paper Data)", x = "SBS4 survival cutoff value", y = "Frequency") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = font_size),
    axis.title = element_text(size = font_size)
  )

# SBS13 TMB
sbs13_tmb_bs_distribution <- ggplot(sbs13_tmb_df, aes(x = cutoff)) +
  geom_histogram(fill = "#7570b3", color = "black", bins = 30) +
  geom_vline(aes(xintercept = SBS13_paper_optimal_TMB_cutoff_for_paper_non_bs$estimate, color = "Non-Bootstrapped Cutoff"), linewidth = 1) +
  geom_vline(aes(xintercept = SBS13_TMB_optimal_cutoff_paper_mean, color = "Bootstrap Mean"), size = 1) +
  scale_color_manual(name = NULL, values = c("Non-Bootstrapped Cutoff" = "red", "Bootstrap Mean" = "green")) +
  scale_x_continuous(limits = c(sbs13_min, sbs13_max)) +
  labs(title = "Bootstrap Distribution of Optimal SBS13 TMB Cutoff (Paper Data)", x = "SBS13 TMB cutoff value", y = "Frequency") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = font_size),
    axis.title = element_text(size = font_size)
  )

# SBS13 Survival
sbs13_surv_bs_distribution <- ggplot(sbs13_surv_df, aes(x = cutoff)) +
  geom_histogram(fill = "#7570b3", color = "black", bins = 30) +
  geom_vline(aes(xintercept = SBS13_paper_optimal_survival_cutoff_for_paper_non_bs$estimate, color = "Non-Bootstrapped Cutoff"), linewidth = 1) +
  geom_vline(aes(xintercept = SBS13_survival_optimal_cutoff_paper_mean, color = "Bootstrap Mean"), size = 1) +
  scale_color_manual(name = NULL, values = c("Non-Bootstrapped Cutoff" = "red", "Bootstrap Mean" = "green")) +
  scale_x_continuous(limits = c(sbs13_min, sbs13_max)) +
  labs(title = "Bootstrap Distribution of Optimal SBS13 Survival Cutoff (Paper Data)", x = "SBS13 survival cutoff value", y = "Frequency") +
  theme_minimal() +
  theme(
    legend.position = "top",
    legend.text = element_text(size = font_size),
    axis.title = element_text(size = font_size)
  )


paper_signatures <- as.data.table(paper_signatures)
paper_signatures_merged <- merge(paper_signatures, biological_weights_paper[, .(Unique_Patient_Identifier, total_snvs)], by = "Unique_Patient_Identifier")
data_paper_sig <- paper_signatures_merged[, .(total_snvs, SBS4, SBS5, SBS13)]
data_paper_sig <- data_paper_sig[complete.cases(data_paper_sig),]

SBS4_0.25_cutoff_paper_signatures <- plot_TMB_binary_boxplot(data_paper_sig, "SBS4", 0.25, "(Binary for Signature Presence: 0.25)")
SBS13_0.25_cutoff_paper_signatures <- plot_TMB_binary_boxplot(data_paper_sig, "SBS13", 0.25, "(Binary for Signature Presence: 0.25)")

# Plot a linear regression of log(TMB) vs Mutations Attributed to SBS4 and SBS13 (paper data)

paper_TMB_linreg_SBS4 <- plot_TMB_linear_regression(biological_weights_paper, "paper", "SBS4")
paper_TMB_linreg_SBS13 <- plot_TMB_linear_regression(biological_weights_paper, "paper", "SBS13")
paper_TMB_linreg_SBS5 <- plot_TMB_linear_regression(biological_weights_paper, "paper", "SBS5")

paper_TMB_linreg_SBS24 <- plot_TMB_linear_regression(biological_weights_paper, "paper", "SBS24")
paper_TMB_linreg_SBS87 <- plot_TMB_linear_regression(biological_weights_paper, "paper", "SBS87")



# Signature Absence (<0.25 or <paper_bs_optimal_cutoff) vs Presence (>= 0.25 or >=paper_bs_optimal_cutoff)

# SBS5_results$paper$binary$plot_0.25
SBS4_survival_0.25 <- SBS4_results$paper$binary$plot_0.25
SBS4_survival_optimal <- SBS4_results$paper$binary$plot_optimal
SBS4_survival_continuous <- SBS4_results$paper$continuous$plot_continuous

SBS13_survival_0.25 <- analyze_signature(signature_name = "SBS13", dataset = "paper", biological_weight_table = biological_weights_paper, survival_data = survival_data_paper, p_value = FALSE)
SBS13_survival_0.25 <- SBS13_survival_0.25$plot_0.25
SBS13_survival_optimal <- SBS13_results$paper$binary$plot_optimal
SBS13_survival_continuous <- SBS13_results$paper$continuous$plot_continuous

SBS4_0.25_cutoff_survival_paper_signatures <- sig_result_SBS4$plot_0.25 # Li et al mutational activities
SBS13_0.25_cutoff_survival_paper_signatures <- sig_result_SBS13$plot_0.25 # Li et al mutational activities



# Mutational Signature Analysis

avg_effect_weights_paper <- signature_attribution_paper$effect_shares$average_effect_shares
diverging_plot_sig_comparison <- create_diverging_signature_plot_with_cancer_effect(Li_data = paper_signatures, ref_based_data = biological_weights_paper, effect_data = avg_effect_weights_paper)



