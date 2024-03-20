###Create Figures

#Figure 1 - plant - general microbial trends
fig1 <- ggarrange(plant.box, pcoa.rel.plot, stack.phylum, labels = c("A","B", "C"), common.legend = F, nrow=2, ncol=2)
fig1
ggsave(fig1, filename = here("output", "figs", "fig1.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 14, units = "in")

#Figure 2 - functional traits relative abundance
fig2 <- ggarrange(ACC.rel.box, IAA.rel.box, osmo.rel.box, EPS.rel.box,
                cyto.rel.box, antiox.rel.box, label.x = c(0.08,0.09,0.040,0.09,0.055,0.035), label.y = 0.98, labels = c("ACC", "IAA", "Osmolytes","EPS", "Cytokinines", "Antioxidants"),
                common.legend = T, nrow=3, ncol=2, legend = "bottom")
fig2
ggsave(fig2, filename = here("output", "figs", "fig2.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 14, units = "in")

#Figure 3 - MAGs affected by actual or historical water stress
fig3.1 <- ggarrange(bar.mag.perc.soil, bar.gene.mag.soil, bar.soil.overlap, labels = c("A","B","C"),
                    common.legend = T, nrow=3, ncol=1, legend = "bottom")
#fig3.2 <- ggarrange(bar.mag.perc.SWHC, bar.gene.mag.SWHC, bar.SWHC.overlap, labels = c("B", "D", "F"),
#                    common.legend = T, nrow=3, ncol=1, legend = "bottom")
#fig3 <- ggarrange(fig3.1,fig3.2, nrow = 1, ncol = 2)

fig3.1
ggsave(fig3.1, filename = here("output", "figs", "fig3.tiff"), compression = "lzw", dpi = 600, device = "tiff", height = 14, width = 14, units = "in")
