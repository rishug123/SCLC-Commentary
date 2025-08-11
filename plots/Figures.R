# setwd("plots") # the location the r script will save

# Fig 1
Fig1a <- diverging_plot_sig_comparison
Fig1b <- SBS4_0.25 +  theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2)) + ggtitle(NULL)
Fig1d <- SBS4_optimal + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2))
Fig1c <- SBS13_0.25 + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2))
Fig1e <- SBS13_optimal + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2))
Fig1f <- paper_TMB_linreg_SBS4 + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2)) 
Fig1g <- paper_TMB_linreg_SBS13 + ggtitle(NULL) + theme(legend.position = "none", legend.margin = margin(0,0,0,0),  legend.spacing.x = unit(1, "mm"), legend.spacing.y = unit(1, "mm")) + theme(aspect.ratio = 0.5, plot.margin = margin(0, 10, 0, 2)) 

Fig1_p1 <- plot_grid(Fig1a, nrow = 1, ncol = 1, labels = "a")
Fig1_p2 <- plot_grid(Fig1b, Fig1c, Fig1d, Fig1e, Fig1f, Fig1g, nrow=3, ncol=2, labels = c("b", "c", "d", "e", "f", "g"))

Fig1 <- plot_grid (Fig1_p1, Fig1_p2, nrow=2, ncol = 1, rel_heights = c(1,1.5))
Fig1 <- Fig1 + theme(plot.background = element_rect(fill = "white", color = "white"))
ggsave("figure_1.png", plot = Fig1, width = 7, height = 10)

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

Fig2e <- SBS4_survival_continuous + theme(
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(1, "mm"),
  legend.spacing.y = unit(1, "mm"),
  aspect.ratio = 1, plot.margin = margin(0, 10, 0, 10)
)
Fig2f <- SBS13_survival_continuous + theme(
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(1, "mm"),
  legend.spacing.y = unit(1, "mm"),
  aspect.ratio = 1, plot.margin = margin(0, 10, 0, 10)
)

panel_a <- plot_grid(Fig2a$plot, Fig2a$table, ncol = 1, align = "v", rel_heights = c(3, 1), labels = c("a", ""),  label_y = 0.98)
panel_b <- plot_grid(Fig2b$plot, Fig2b$table, ncol = 1, align = "v", rel_heights = c(3, 1), labels = c("b", ""),  label_y = 0.98)
panel_c <- plot_grid(Fig2c$plot, Fig2c$table, ncol = 1, align = "v", rel_heights = c(3, 1), labels = c("c", ""),  label_y = 0.98)
panel_d <- plot_grid(Fig2d$plot, Fig2d$table, ncol = 1, align = "v", rel_heights = c(3, 1), labels = c("d", ""),  label_y = 0.98)

# Combine all panels into a 2x2 grid with tags A-D
Fig2 <- plot_grid(panel_a, panel_b, panel_c, panel_d, Fig2e, Fig2f, nrow=3, ncol=2, labels = c("", "", "", "", "e", "f"), rel_heights = c(1, 1, 0.8)) + theme(plot.background = element_rect(fill = "white", color = "white"),  label_y = 2, vjust = -0.5)
ggsave("figure_2.png", plot = Fig2, width = 7, height = 12)


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

# Assign the external survival plots properly

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

sfig3c <- SBS4_results$external$continuous$plot_continuous + theme(
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(1, "mm"),
  legend.spacing.y = unit(1, "mm"),
  aspect.ratio = 1, plot.margin = margin(0, 10, 0, 10)
)
sfig3d <- SBS13_results$external$continuous$plot_continuous + theme(
  legend.margin = margin(0, 0, 0, 0),
  legend.spacing.x = unit(1, "mm"),
  legend.spacing.y = unit(1, "mm"),
  aspect.ratio = 1, plot.margin = margin(0, 10, 0, 10)
)

# Combine into 2x2 grid with tags
SuppFig3 <- plot_grid(sfig3a, sfig3b, sfig3c, sfig3d, 
                      nrow = 2, ncol = 2, labels = c("a", "b", "c", "d"), rel_heights = c(1,1)) +
  theme(plot.background = element_rect(fill = "white", color = "white"))

# Save figure
ggsave("supp_fig_3.png", plot = SuppFig3, width = 7, height = 7)


