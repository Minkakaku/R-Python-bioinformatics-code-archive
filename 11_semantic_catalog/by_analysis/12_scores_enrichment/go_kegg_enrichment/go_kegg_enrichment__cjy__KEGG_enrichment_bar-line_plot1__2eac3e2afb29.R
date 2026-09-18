library(readxl)
library(dplyr)
library(ggplot2)

kegg.enrich <- read_excel("dow.xlsx")

kegg.enrich <- kegg.enrich %>%
  mutate(
    pvalue = as.numeric(pvalue),
    Count = as.numeric(Count),
    logP = -log10(pvalue),
    Count_scaled = Count / max(Count, na.rm = TRUE) * max(logP, na.rm = TRUE)
  )

p <- ggplot(
  kegg.enrich,
  aes(x = logP, y = reorder(Description, pvalue))
) +
  geom_col(fill = "pink", width = 0.7) +
  geom_point(aes(x = Count_scaled), color = "black", size = 2) +
  geom_path(
    aes(x = Count_scaled, y = reorder(Description, pvalue), group = 1),
    color = "grey",
    linewidth = 1
  ) +
  # scale_x_continuous(
  #   expand = c(0, 0),
  #   name = "-log10(pvalue)",
  #   sec.axis = sec_axis(
  #     ~ . / max(kegg.enrich$logP, na.rm = TRUE) * max(kegg.enrich$Count, na.rm = TRUE),
  #     name = "Count"
  #   )
  # ) +
  scale_x_continuous(
    expand = expansion(mult = c(0, 0.08)),
    name = "-log10(pvalue)",
    sec.axis = sec_axis(
      ~ . / max(kegg.enrich$logP, na.rm = TRUE) * max(kegg.enrich$Count, na.rm = TRUE),
      name = "Count"
    )
  ) +
  theme_classic(base_size = 14) +
  labs(y = "", title = "KEGG enrichment barplot") +
  theme(
    legend.title = element_text(hjust = 0.5)
  )

p

ggsave("KEGG_enrichment_barline.pdf", p, width = 8, height = 6)
ggsave("KEGG_enrichment_barline.png", p, width = 8, height = 6, dpi = 300)