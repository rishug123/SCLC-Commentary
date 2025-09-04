setwd("plots") # the location the r script will save

# Fig 1
Fig1a_1 <- diverging_plot_sig_comparison$TMB
Fig1a_2 <- diverging_plot_sig_comparison$label + theme(legend.margin = margin(0, 0, 0, 0))
Fig1a_3 <- diverging_plot_sig_comparison$effect


Fig1b <- SBS4_0.25 +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 10, 0, 2),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )
Fig1c <- SBS13_0.25 +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 10, 0, 2),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )


# SBS4 survival at 0.25 cutoff
Fig1d <- SBS4_survival_0.25
Fig1d$plot <- Fig1d$plot +
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
Fig1e <- SBS13_survival_0.25
Fig1e$plot <- Fig1e$plot +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "top",
    legend.margin = margin(0,0,0,0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )



panel_d <- plot_grid(Fig1d$plot, Fig1d$table, ncol = 1, rel_heights = c(3, 1), labels = c("", ""),  label_y = 0.98)
panel_e <- plot_grid(Fig1e$plot, Fig1e$table, ncol = 1,  rel_heights = c(3, 1), labels = c("", ""),  label_y = 0.98)


Fig1_p1 <- plot_grid(Fig1a_1, NULL ,Fig1a_3, nrow = 3, ncol = 1, labels = c("a","", "b"), rel_heights = c(12, 1, 7),  label_y = 1.04)
Fig1_top <- plot_grid(Fig1_p1, Fig1a_2, nrow=1, ncol = 2, rel_widths = c(6,4))
Fig1_bottom <- plot_grid(Fig1b, Fig1c, panel_d, panel_e, nrow=1, ncol=4, labels = c("c", "d", "e", "f"),  label_y = 0.97) + theme(plot.background = element_rect(fill = "white", color = "white"))

Fig1 <- plot_grid(Fig1_top, Fig1_bottom, nrow=2, rel_heights = c(1,1)) + theme(plot.background = element_rect(fill = "white", color = "white"))
ggsave("figure_1.png", plot = Fig1, width = 15, height = 8)

# Fig 2

Fig2a <- SBS4_optimal +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 10, 0, 2),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )
Fig2b <- SBS13_optimal +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 10, 0, 2),
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

Fig2c <- plot_grid(Fig2c$plot, Fig2c$table, ncol = 1,  rel_heights = c(3, 1), labels = c("", ""),  label_y = 0.98)
Fig2d <- plot_grid(Fig2d$plot, Fig2d$table, ncol = 1,  rel_heights = c(3, 1), labels = c("", ""),  label_y = 0.98)


Fig2e <- paper_TMB_linreg_SBS4 +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

Fig2f <- paper_TMB_linreg_SBS13 +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )


# Combine into 2x2 grid with tags
Fig2 <- plot_grid(Fig2a, Fig2b, Fig2c, Fig2d, Fig2e, Fig2f, 
                      nrow = 3, ncol = 2, labels = c("a", "b", "c", "d", "e", "f"), rel_heights = c(1,1)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# Save figure
ggsave("figure_2.png", plot = Fig2, width = 7, height = 11)


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

panel_s1c <- plot_grid(sfig1c$plot, sfig1c$table, ncol = 1, align = "v", rel_heights = c(9, 3), labels = c("c", ""),  label_y = 0.96)
panel_s1d <- plot_grid(sfig1d$plot, sfig1d$table, ncol = 1, align = "v", rel_heights = c(9, 3), labels = c("d", ""),  label_y = 0.96)

# Combine into 2x2 grid with tags
SuppFig1 <- plot_grid(sfig1a, sfig1b, panel_s1c, panel_s1d, 
                      nrow = 2, ncol = 2, labels = c("a", "b"), rel_heights = c(2,3)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# Save the figure
ggsave("supp_fig_1.png", plot = SuppFig1, width = 7, height = 8.5)

# supplemental figure 2

sfig2c <- sbs4_tmb_bs_distribution +
  ggtitle(NULL)
sfig2d <- sbs13_tmb_bs_distribution +
  ggtitle(NULL)
sfig2a <- sbs4_surv_bs_distribution +
  ggtitle(NULL)
sfig2b <- sbs13_surv_bs_distribution +
  ggtitle(NULL)

SuppFig2 <- plot_grid(sfig2a, sfig2b, 
                      nrow = 1, ncol = 2, labels = c("a", "b"),  label_y = 0.98)+
  theme(plot.background = element_rect(fill = "white", color = "white"))

ggsave("supp_fig_2.png", plot = SuppFig2, width = 8, height = 3.5)

# supplemental figure 3

sfig3a <- external_TMB_linreg_SBS4
sfig3b <- external_TMB_linreg_SBS13

# Apply consistent theme to each plot in Supplemental Figure 1
sfig3a <- sfig3a +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig3b <- sfig3b +
  ggtitle(NULL) +
  labs(color = NULL) +
  theme(
    legend.position = "none",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

sfig3c <- SBS4_results$external$binary$plot_optimal
sfig3d <- SBS13_results$external$binary$plot_optimal

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

sfig3d$plot <- sfig3d$plot +
  ggtitle(NULL) +
  labs(color = NULL)+
  theme(
    legend.position = "top",
    legend.margin = margin(0, 0, 0, 0),
    legend.spacing.x = unit(1, "mm"),
    legend.spacing.y = unit(1, "mm"),
    aspect.ratio = 1
  )

panel_s3c <- plot_grid(sfig3c$plot, sfig3c$table, ncol = 1, align = "v", rel_heights = c(9, 3), labels = c("c", ""),  label_y = 0.96)
panel_s3d <- plot_grid(sfig3d$plot, sfig3d$table, ncol = 1, align = "v", rel_heights = c(9, 3), labels = c("d", ""),  label_y = 0.96)



# Combine into 2x2 grid with tags
SuppFig3 <- plot_grid(sfig3a, sfig3b, panel_s3c, panel_s3d, 
                      nrow = 2, ncol = 2, labels = c("a", "b", "", ""), rel_heights = c(2,3)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# Save figure
ggsave("supp_fig_3.png", plot = SuppFig3, width = 7, height = 7)


