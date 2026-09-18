# ============================================================
# GO BP / CC / MF GeneRatio calculation + facet dot plot
# Based on existing enrichment results
# ============================================================


library(data.table)
library(dplyr)
library(openxlsx)
library(ggplot2)
library(stringr)



# ============================================================
# 1. 读取 GO annotation Table3
# ============================================================

bp_annotation <- fread(
  file.choose(),
  check.names = FALSE
)

cc_annotation <- fread(
  file.choose(),
  check.names = FALSE
)

mf_annotation <- fread(
  file.choose(),
  check.names = FALSE
)



# ============================================================
# 2. 读取 GO enrichment Table4
# ============================================================

bp_enrich <- fread(
  file.choose(),
  check.names = FALSE
)

cc_enrich <- fread(
  file.choose(),
  check.names = FALSE
)

mf_enrich <- fread(
  file.choose(),
  check.names = FALSE
)



# ============================================================
# 3. 标准化列名
# ============================================================

clean_col <- function(df){
  
  colnames(df) <- gsub("\\.", " ", colnames(df))
  
  idx <- grep(
    "accession",
    colnames(df),
    ignore.case = TRUE
  )
  
  if(length(idx) > 0){
    colnames(df)[idx[1]] <- "Protein accession"
  }
  
  return(df)
}



bp_annotation <- clean_col(bp_annotation)
cc_annotation <- clean_col(cc_annotation)
mf_annotation <- clean_col(mf_annotation)


bp_enrich <- clean_col(bp_enrich)
cc_enrich <- clean_col(cc_enrich)
mf_enrich <- clean_col(mf_enrich)



# ============================================================
# 4. 计算各GO类别 GeneRatio
# ============================================================

calc_gene_ratio <- function(annotation,
                            enrich,
                            class_name){
  
  
  # 注释成功的蛋白数量
  n_background <- annotation %>%
    select(`Protein accession`) %>%
    distinct() %>%
    nrow()
  
  
  message(
    class_name,
    " annotated proteins = ",
    n_background
  )
  
  
  enrich %>%
    mutate(
      mapping = as.numeric(mapping),
      
      GeneRatio = mapping / n_background,
      
      GO_class = class_name
    )
  
}



bp_result <- calc_gene_ratio(
  bp_annotation,
  bp_enrich,
  "GO-BP"
)


cc_result <- calc_gene_ratio(
  cc_annotation,
  cc_enrich,
  "GO-CC"
)


mf_result <- calc_gene_ratio(
  mf_annotation,
  mf_enrich,
  "GO-MF"
)



# ============================================================
# 5. 每类按p value筛选Top10
# ============================================================

bp_top <- bp_result %>%
  arrange(`p value`) %>%
  slice_head(n = 10)


cc_top <- cc_result %>%
  arrange(`p value`) %>%
  slice_head(n = 10)


mf_top <- mf_result %>%
  arrange(`p value`) %>%
  slice_head(n = 10)



# 合并
go_plot_data <- bind_rows(
  bp_top,
  cc_top,
  mf_top
)



# ============================================================
# 6. 输出GeneRatio结果表
# ============================================================

write.xlsx(
  go_plot_data,
  "GO_GeneRatio_Top10.xlsx"
)



# ============================================================
# 7. 设置Term排序
# ============================================================

go_plot_data$Terms <- factor(
  go_plot_data$Terms,
  levels =
    go_plot_data %>%
    arrange(GeneRatio) %>%
    pull(Terms)
)



# ============================================================
# 8. GO facet dot plot
# ============================================================

p <- ggplot(
  go_plot_data,
  aes(
    x = GeneRatio,
    y = Terms
  )
) +
  
  geom_point(
    aes(
      size = mapping,
      color = `p value`
    )
  ) +
  
  facet_grid(
    GO_class ~ .,
    scales = "free_y",
    space = "free_y"
  ) +
  
  theme_bw(
    base_size = 14
  ) +
  
  labs(
    x = "Gene Ratio",
    y = "",
    title = "",
    color = "p value",
    size = "Protein number"
  ) +
  
  scale_y_discrete(
    labels = function(x){
      stringr::str_wrap(
        x,
        width = 60
      )
    }
  ) +
  
  scale_color_gradientn(
    colors = c(
      "red",
      "orange",
      "yellow",
      "green3",
      "skyblue",
      "darkblue"
    )
  ) +
  
  scale_size_continuous(
    range = c(3,8)
  ) +
  
  scale_x_continuous(
    expand = expansion(
      mult = c(0.05,0.15)
    )
  ) +
  
  theme(
    
    axis.text = element_text(
      size = 12.5
    ),
    
    panel.grid.major = element_blank(),
    
    panel.grid.minor = element_blank(),
    
    strip.background = element_blank(),
    
    strip.text.y = element_text(
      size = 12,
      face = "bold"
    ),
    
    legend.title = element_text(
      face = "bold"
    )
  )



# ============================================================
# 9. 保存图片
# ============================================================

ggsave(
  "GO_Facet_GeneRatio_dotplot.pdf",
  p,
  width = 10,
  height = 10
)


ggsave(
  "GO_Facet_GeneRatio_dotplot.png",
  p,
  width = 10,
  height = 10,
  dpi = 300
)



# 显示图片
print(p)