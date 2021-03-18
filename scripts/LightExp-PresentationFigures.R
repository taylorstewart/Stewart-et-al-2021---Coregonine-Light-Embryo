source("scripts/LightExp-Embryo-LHT.R")

## Combine all figures
plot.survival.all <- grid.arrange(
  arrangeGrob(
    arrangeGrob(textGrob(""),
                get_legend(plot.survival + theme(legend.key.size = unit(0.75, 'cm'))),
                nrow = 1,
                widths = c(0.09, 1)),
    arrangeGrob(textGrob(""),
                get_legend(plot.survival.stand + theme(legend.key.size = unit(0.75, 'cm'))),
                nrow = 1,
                widths = c(0.09, 1)),
    ncol = 2,
    widths = c(0.7, 0.7)
  ),
  arrangeGrob(
    arrangeGrob(plot.survival + theme(legend.position = "none", axis.title.x = element_blank()),
                bottom = textGrob("Population", x = 0.6, gp = gpar(cex = 2, fontfamily = "Arial"))),
    arrangeGrob(plot.survival.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                bottom = textGrob("Population", x = 0.6, gp = gpar(cex = 2, fontfamily = "Arial"))), 
    ncol = 2,
    widths = c(0.7, 0.7)
  ),
  heights = c(0.065, 1.0)
)

ggsave("figures/presentation/survival.png", plot = plot.survival.all, width = 10, height = 5.5, dpi = 200)


## Combine all figures
plot.inc.all <- grid.arrange(
  arrangeGrob(
    arrangeGrob(textGrob(""),
                get_legend(plot.dpf),
                nrow = 1,
                widths = c(0.09, 1)),
    arrangeGrob(textGrob(""),
                get_legend(plot.dpf.stand + theme(legend.key.size = unit(1, 'cm'))),
                nrow = 1,
                widths = c(0.09, 1)),
    ncol = 2,
    widths = c(0.7, 0.7)
  ),
  arrangeGrob(
    arrangeGrob(
      plot.dpf + theme(legend.position = "none", axis.title.x = element_blank()) +
        labs(y = "Mean Days Post-Fertilization"),
      plot.ADD + theme(legend.position = "none", axis.title.x = element_blank()) +
        labs(y = "Mean Accumulated Degree Days"),
      nrow = 2,
      bottom = textGrob("Population", x = 0.61, gp = gpar(cex = 2, fontfamily = "Arial"))
    ),
    arrangeGrob(plot.dpf.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.ADD.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 2,
                bottom = textGrob("Population", x = 0.64, gp = gpar(cex = 2, fontfamily = "Arial")),
                ncol = 1),
    widths = c(0.7, 0.7)
  ),
  heights = c(0.055, 1.0)
)

ggsave("figures/presentation/DPF-ADD.png", plot = plot.inc.all, width = 10, height = 10, dpi = 200)



source("scripts/LightExp-Larvae-MT.R")

## Combine all figures
plot.mt.all <- grid.arrange(
  arrangeGrob(
    arrangeGrob(textGrob(""),
                get_legend(plot.tl),
                nrow = 1,
                widths = c(0.09, 1)),
    arrangeGrob(textGrob(""),
                get_legend(plot.tl.stand + theme(legend.key.size = unit(1, 'cm'))),
                nrow = 1,
                widths = c(0.09, 1)),
    ncol = 2,
    widths = c(0.7, 0.7)
  ),
  arrangeGrob(
    arrangeGrob(
      plot.tl + theme(legend.position = "none", axis.title.x = element_blank()) +
        labs(y = "Mean Length-at-Hatch (mm)"),
      plot.yolk + theme(legend.position = "none", axis.title.x = element_blank()) +
        labs(y = expression("Mean Yolk-sac Volume (mm"^3*")")),
      nrow = 2,
      bottom = textGrob("Population", x = 0.61, gp = gpar(cex = 2, fontfamily = "Arial"))
    ),
    arrangeGrob(plot.tl.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                plot.yolk.stand + theme(legend.position = "none", axis.title.x = element_blank()),
                nrow = 2,
                bottom = textGrob("Population", x = 0.64, gp = gpar(cex = 2, fontfamily = "Arial")),
                ncol = 1),
    widths = c(0.7, 0.7)
  ),
  heights = c(0.055, 1.0)
)


ggsave("figures/presentation/tl-ysv.png", plot = plot.mt.all, width = 10, height = 10, dpi = 200)

