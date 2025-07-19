setwd("plots") # the location the r script will save

# Fig 1
Fig1a <- diverging_plot_sig_comparison
Fig1d <- SBS4_0.25 + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2)) 
Fig1f <- SBS4_optimal + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2))
Fig1e <- SBS13_0.25 + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2))
Fig1g <- SBS13_optimal + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2))
Fig1b <- paper_TMB_linreg_SBS4 + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2)) 
Fig1c <- paper_TMB_linreg_SBS13 + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2)) 

Fig1_p1 <- plot_grid(Fig1a, nrow = 1, ncol = 1, labels = "a")
Fig1_p2 <- plot_grid(Fig1b, Fig1c, Fig1d, Fig1e, Fig1f, Fig1g, nrow=3, ncol=2, labels = c("b", "c", "d", "e", "f", "g"))

Fig1 <- plot_grid (Fig1_p1, Fig1_p2, nrow=2, ncol = 1, rel_heights = c(1,2))
Fig1 <- Fig1 + theme(plot.background = element_rect(fill = "white", color = "white"))
ggsave("figure_1.png", plot = Fig1, width = 7, height = 9)

# Fig 2

# SBS4 survival at 0.25 cutoff
Fig2a <- SBS4_survival_0.25
Fig2a$plot <- Fig2a$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0,0,0,0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

# SBS13 survival at 0.25 cutoff
Fig2b <- SBS13_survival_0.25
Fig2b$plot <- Fig2b$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0,0,0,0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

# SBS4 survival at optimal cutoff
Fig2c <- SBS4_survival_optimal
Fig2c$plot <- Fig2c$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0,0,0,0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

# SBS13 survival at optimal cutoff
Fig2d <- SBS13_survival_optimal
Fig2d$plot <- Fig2d$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0,0,0,0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

panel_a <- plot_grid(Fig2a$plot, Fig2a$table, ncol = 1, align = "v", rel_heights = c(3, 1))
panel_b <- plot_grid(Fig2b$plot, Fig2b$table, ncol = 1, align = "v", rel_heights = c(3, 1))
panel_c <- plot_grid(Fig2c$plot, Fig2c$table, ncol = 1, align = "v", rel_heights = c(3, 1))
panel_d <- plot_grid(Fig2d$plot, Fig2d$table, ncol = 1, align = "v", rel_heights = c(3, 1))

# Combine all panels into a 2x2 grid with tags A-D
Fig2 <- plot_grid(panel_a, panel_b, panel_c, panel_d, nrow=2, ncol=2, labels = c("a", "b", "c", "d")) + theme(plot.background = element_rect(fill = "white", color = "white"))
ggsave("figure_2.png", plot = Fig2, width = 7, height = 8.5)


# Supp fig 1

sfig1a <- SBS4_0.25_cutoff_paper_signatures
sfig1b <- SBS13_0.25_cutoff_paper_signatures
sfig1c <- SBS4_0.25_cutoff_survival_paper_signatures
sfig1d <- SBS13_0.25_cutoff_survival_paper_signatures

# Apply consistent theme to each plot in Supplemental Figure 1
sfig1a <- sfig1a +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig1b <- sfig1b +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig1c$plot <- sfig1c$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig1d$plot <- sfig1d$plot +
  ggtitle(NULL) +
  labs(color = NULL)+
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

# Stack plot and table for each panel

panel_s1c <- plot_grid(sfig1c$plot, sfig1c$table, ncol = 1, align = "v", rel_heights = c(3, 1))
panel_s1d <- plot_grid(sfig1d$plot, sfig1d$table, ncol = 1, align = "v", rel_heights = c(3, 1))

# Combine into 2x2 grid with tags
SuppFig1 <- plot_grid(sfig1a, sfig1b, panel_s1c, panel_s1d, 
                      nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"), rel_heights = c(2,3)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# Save the figure
ggsave("supp_fig_1.png", plot = SuppFig1, width = 7, height = 8.5)

# supplemental figure 2

sfig2a <- sbs4_tmb_bs_distribution +
  ggtitle(NULL)
sfig2b <- sbs13_tmb_bs_distribution +
  ggtitle(NULL)
sfig2c <- sbs4_surv_bs_distribution +
  ggtitle(NULL)
sfig2d <- sbs13_surv_bs_distribution +
  ggtitle(NULL)

SuppFig2 <- plot_grid(sfig2a, sfig2b, sfig2c, sfig2d, 
                      nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"))+
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("supp_fig_2.png", plot = SuppFig2, width = 8, height = 7)

# Assign the external survival plots properly
SBS4_survival_0.25_external <- SBS4_results$external$plot_0.25
SBS13_survival_0.25_external <- SBS13_results$external$plot_0.25
SBS4_survival_optimal_external <- SBS4_results$external$plot_optimal
SBS13_survival_optimal_external <- SBS13_results$external$plot_optimal

# Format plots
sfig3a <- SBS4_survival_0.25_external
sfig3a$plot <- sfig3a$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig3b <- SBS13_survival_0.25_external
sfig3b$plot <- sfig3b$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig3c <- SBS4_survival_optimal_external
sfig3c$plot <- sfig3c$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig3d <- SBS13_survival_optimal_external
sfig3d$plot <- sfig3d$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

# Stack plot and table for each subplot
panel_s3a <- plot_grid(sfig3a$plot, sfig3a$table, ncol = 1, align = "v", rel_heights = c(3, 1))
panel_s3b <- plot_grid(sfig3b$plot, sfig3b$table, ncol = 1, align = "v", rel_heights = c(3, 1))
panel_s3c <- plot_grid(sfig3c$plot, sfig3c$table, ncol = 1, align = "v", rel_heights = c(3, 1))
panel_s3d <- plot_grid(sfig3d$plot, sfig3d$table, ncol = 1, align = "v", rel_heights = c(3, 1))

# Combine into 2x2 grid with panel labels
SuppFig3 <- plot_grid(panel_s3a, panel_s3b, panel_s3c, panel_s3d,
                      nrow = 2, ncol = 2,
                      labels = c("a", "b", "c", "d")) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# Save figure
ggsave("supp_fig_3.png", plot = SuppFig3, width = 7, height = 8.5)


