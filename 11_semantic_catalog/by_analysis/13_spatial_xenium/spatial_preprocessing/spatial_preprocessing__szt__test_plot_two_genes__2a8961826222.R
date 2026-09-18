rownames(all_data) %>% head
# "Xkr4"   "Rp1"    "Sox17"  "Lypla1" "Tcea1"  "Rgs20" 

# HE 
p <- SpatialDimPlot(
  all_data,
  images = 'SYZ_NH4.008um',
  # features = NULL,
  group.by = NULL,
  alpha = 0
) 

p$data %>% head



dat_gene <- FetchData(all_data, vars = c("Mbp",  "mt-Co2"))


dat_gene_full <- 
  GetTissueCoordinates(all_data, scale = "lowres", image = 'SYZ_NH4.008um') %>%
  left_join(
    rownames_to_column(dat_gene, "cell"),
    by = "cell"
  )  

dat_gene_full_ls <-
  dat_gene_full %>%
  mutate(
    Mbp2 = as.numeric(scale(Mbp)),
    `mt-Co2-2` = as.numeric(scale(`mt-Co2`))
  ) %>%
  mutate(group = ifelse(Mbp2 > `mt-Co2-2`, 'Mbp', 'mt-Co2')) %>%
  group_by(group) %>%
  nest() %>%
  deframe()  

# rownames(all_data)[c(17651, 19050, 19058, 19052, 19051, 19049)]


ggplot() +
  geom_point(
    data = dat_gene_full_ls$Mbp %>% mutate(x = max(x) - x + min(x)),
    aes(x = y, y = x, color = Mbp), #, height = 0.2, width = 0.2),
    size = 0.1, 
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(colours = c("white", "gray90", "blue")) +
  ggnewscale::new_scale_color() +
  geom_point(
    data = dat_gene_full_ls$`mt-Co2` %>% mutate(x = max(x) - x + min(x)),
    aes(x = y, y = x, color = `mt-Co2`),
    size = 0.1, 
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  scale_color_gradientn(colours = c("white", "gray90", "red")) +
  theme_void() +
  coord_fixed()



dat <-
  data.frame(
    x = 1:5, 
    y = 1:5
  )


dat %>%
  ggplot() +
  geom_tile(aes(x = x, y = y, height = 0.2, width = 0.4)) +
  geom_tile(aes(x = x, y = y, height = 0.1, width = 0.2), fill = "red",   position = position_nudge(y = -0.05)) +
  geom_tile(aes(x = x, y = y, height = 0.1, width = 0.2), fill = "blue",  position = position_nudge(y =  0.05))


# new 
height_bin <- 0.25 # 调整height大小，保证不遮挡即可
width_bin  <- 0.5 # 调整height大小，保证不遮挡即可
ggplot() +
  geom_tile(
    data = dat_gene_full %>% mutate(x = max(x) - x + min(x)),
    aes(x = y, y = x, fill = Mbp, height = height_bin, width = width_bin),
    position = position_nudge(y = - height_bin /2),
    size = 0.1, 
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  scale_fill_gradientn(colours = c("white", "blue")) +
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = dat_gene_full %>% mutate(x = max(x) - x + min(x)),
    aes(x = y, y = x, fill = `mt-Co2`, height = height_bin, width = width_bin),
    position = position_nudge(y = height_bin/2),
    size = 0.1, 
    alpha = 0.6,
    inherit.aes = FALSE
  ) +
  scale_fill_gradientn(colours = c("white", "red")) +
  theme_void() +
  coord_fixed()
ggsave('test_two_gene.pdf', height = 10, width = 10)
  