library(ggplot2)

setwd('C:/Users/ztt/Desktop/niches')

data = read.table('proportion.txt', sep='\t', header = T)

nuclei_niche_color_list <- c(
  `C0` = "#1036b5",
  `C1` = "#6f774b",
  `C2` = "#d10040",
  `C3` = "#84cd5d",
  `C4` = "#8E76DB",
  `C5` = "#71c8a5",
  `C6` = "#ffa26d",
  `C7` = "#924373",
  `C8` = "#a6513c",
  `C9` = "black",
  `C10` = "#72c7ff",
  `C11` = "#cabd00"
)

ggplot(data, aes(x=celltype, y=niches, size=Proportion, fill = as.factor(niches))) +
  geom_point(shape = 21, color = "black") +
  # guides(fill = "none", size = guide_legend(override.aes = list(fill = "red")))
  # guides(fill = "none") +
  labs(title = 'Cellular niches') +
  theme_bw() +
  scale_fill_manual(values = nuclei_niche_color_list, guide = "none") +
  scale_color_manual(values = c("black", "grey80", rep("black", 10)), guide = "none") +
  scale_size_continuous(range = c(0, 10), breaks = seq(0.2, 1, 0.2), limits = c(0, 1)) +
  theme(axis.text = element_text(color = "black", face='bold', size=12),
        axis.text.x = element_text(angle = 90, hjust = 1, vjust = 0.5),
        panel.grid.major.x = element_blank(),
        axis.title = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.y = element_line(color = "grey90", linewidth = 0.3),
        strip.background = element_blank(), strip.text = element_blank(),
        plot.title = element_text(color = "black", face='bold', size=14, hjust = 0.5),
        plot.margin = margin(5.5, 5.5, 5.5, 2, "pt"),
        legend.key = element_blank(),
        legend.box.margin = margin(0, 0, 0, 0, "pt"),
        legend.background = element_blank(),
        legend.direction = "horizontal",
        legend.text = element_text(margin = margin(0, 0, 0, 0, "pt")),
        legend.spacing = unit(0, "pt"),
        legend.position = "bottom", legend.title.align = 0.5,
        legend.key.width = unit(3.5, "mm"),
        legend.key.height = unit(2, "mm")) +
  guides(size = guide_legend(title.position = "top", label.position = "bottom")) +
  scale_alpha_manual(values = c(0, 1), guide = "none")


