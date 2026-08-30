# StrainPhlAn t__SGB7222 Lentilactobacillus parakefiri Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/tree_and_phylogenetic_distance
# 17/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkarian' = '#4FAFA8', 'Karachay-Cherkess' = '#E8403D')


# Load the data -----------------------------------------------------------

lparakefiri_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7222_Lentilactobacillus_parakefiri_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
lparakefiri_mat_df <- lparakefiri_data |>
  bind_rows(
    lparakefiri_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
lparakefiri_all_ids <- sort(unique(c(lparakefiri_mat_df$strain_id_1, lparakefiri_mat_df$strain_id_2)))

lparakefiri_diag_df <- tibble(
  strain_id_1 = lparakefiri_all_ids,
  strain_id_2 = lparakefiri_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
lparakefiri_mat_df <- bind_rows(lparakefiri_diag_df, lparakefiri_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)


# Heatmap -----------------------------------------------------------------

# Create region annotation table (one row per strain)
lparakefiri_region_df <- tibble(strain_id = lparakefiri_all_ids) |>
  mutate(
    collection_site = if_else(
      strain_id %in% kab_samples,
      'Kabardino-Balkarian',
      'Karachay-Cherkess'
    ),
    collection_site = factor(
      collection_site,
      levels = c('Kabardino-Balkarian', 'Karachay-Cherkess')
    )
  )

# Keep only sites that actually exist in this dataset (prevents KC showing in legend when absent)
lparakefiri_sites_present <- lparakefiri_region_df |>
  distinct(collection_site) |>
  pull(collection_site) |>
  as.character()

# Convert heatmap coordinates to numeric indices (so we can place annotations cleanly)
lparakefiri_index_df <- tibble(
  strain_id = lparakefiri_all_ids,
  x = seq_along(lparakefiri_all_ids)
)

lparakefiri_index_df_y <- tibble(
  strain_id = rev(lparakefiri_all_ids),
  y = seq_along(rev(lparakefiri_all_ids))
)

lparakefiri_plot_df <- lparakefiri_mat_df |>
  left_join(lparakefiri_index_df, by = c('strain_id_1' = 'strain_id')) |>
  left_join(lparakefiri_index_df_y, by = c('strain_id_2' = 'strain_id')) |>
  rename(x = x, y = y)

# Annotation positions (outside tile area, but NOT stretching axes)
lparakefiri_n <- length(lparakefiri_all_ids)

lparakefiri_hm_left_x   <- 0.25  # left of first tile center (tile centers start at 1)
lparakefiri_hm_bottom_y <- 0.25  # below the last tile row


# Plot the data -----------------------------------------------------------

lparakefiri_hm <- ggplot(lparakefiri_plot_df, aes(x = x, y = y, fill = phylogenetic_distance)) +
  geom_tile() +
  
  # Add distance values inside each tile
  geom_text(
    aes(label = sprintf('%.2f', phylogenetic_distance)),
    color = 'white',
    size  = 3.5
  ) +
  
  # LEFT annotation (Caucasus Region) aligned with rows
  geom_point(
    data = lparakefiri_region_df |>
      left_join(lparakefiri_index_df_y |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = lparakefiri_hm_left_x,
      y = y,
      shape = collection_site,
      color = collection_site
    ),
    inherit.aes = FALSE,
    size = 3,
    stroke = 0
  ) +
  
  # BOTTOM annotation (Caucasus Region) aligned with columns
  geom_point(
    data = lparakefiri_region_df |>
      left_join(lparakefiri_index_df |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = x,
      y = lparakefiri_hm_bottom_y,
      shape = collection_site,
      color = collection_site
    ),
    inherit.aes = FALSE,
    size = 3,
    stroke = 0
  ) +
  
  # Shapes and colors for Caucasus Region
  scale_shape_manual(
    name   = 'Caucasus\nRegion',
    values = c('Kabardino-Balkarian' = 17, 'Karachay-Cherkess' = 15),
    breaks = lparakefiri_sites_present,
    drop   = TRUE
  ) +
  scale_color_manual(
    name   = 'Caucasus\nRegion',
    values = site_cols,
    breaks = lparakefiri_sites_present,
    drop   = TRUE
  ) +
  
  # Heatmap fill scale
  scale_fill_gradient(
    name = 'Phylogenetic\nDistance',
    low  = '#081d3a',
    high = '#56b1f7',
    limits = c(0, max(lparakefiri_plot_df$phylogenetic_distance, na.rm = TRUE)),
    na.value = 'white'
  ) +
  
  # Axis labels using numeric breaks -> strain IDs
  scale_x_continuous(
    breaks = seq_len(lparakefiri_n),
    labels = lparakefiri_all_ids,
    expand = c(0, 0),
    limits = c(0.0, lparakefiri_n + 0.5)
  ) +
  scale_y_continuous(
    breaks = seq_len(lparakefiri_n),
    labels = rev(lparakefiri_all_ids),
    expand = c(0, 0),
    limits = c(lparakefiri_n + 0.5, 0.0)
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    axis.title   = element_blank(),
    axis.text = element_text(size = 10),
    axis.text.x  = element_text(angle = 90, vjust = 0.5, hjust = 1),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text  = element_text(size = 12)
  ) +
  
  # Legend order: Caucasus Region ABOVE Phylogenetic Distance
  # IMPORTANT: override.aes must match the number of legend keys (sites present), or it errors
  guides(
    shape = 'none',
    color = 'none',
    fill  = guide_colorbar(order = 1)
  )

# Check the plot
lparakefiri_hm


# Making the Script Unique for L. parakefiri ------------------------------

# Load the necessary packages.

library(ape)
library(Biostrings)
library(paletteer)
library(patchwork)
library(ggnewscale)


# Load the data -----------------------------------------------------------

lparakefiri_tree_file <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7222_Lentilactobacillus_parakefiri/RAxML_bestTree.t__SGB7222.StrainPhlAn4.tre'
lparakefiri_msa_file  <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7222_Lentilactobacillus_parakefiri/t__SGB7222.StrainPhlAn4_concatenated.aln'

# Read tree + MSA
lparakefiri_tr  <- read.tree(lparakefiri_tree_file)
lparakefiri_msa <- readDNAStringSet(lparakefiri_msa_file, format = 'fasta') # width = 64

# Match/reorder
lparakefiri_keep <- intersect(lparakefiri_tr$tip.label, names(lparakefiri_msa))
if (length(lparakefiri_keep) == 0) stop('No overlap between tree tip labels and MSA names.')

lparakefiri_tr  <- keep.tip(lparakefiri_tr, lparakefiri_keep)
lparakefiri_msa <- lparakefiri_msa[lparakefiri_tr$tip.label]


# Define clusters (k <= 10 for Classic_10) --------------------------------

lparakefiri_k <- 2
if (lparakefiri_k > 10) stop('k must be <= 10 when using ggthemes::Classic_10.')

lparakefiri_d  <- cophenetic.phylo(lparakefiri_tr)
lparakefiri_hc <- hclust(as.dist(lparakefiri_d), method = 'average')
lparakefiri_cl <- cutree(lparakefiri_hc, k = lparakefiri_k)
lparakefiri_cl <- setNames(as.integer(lparakefiri_cl), names(lparakefiri_cl))

lparakefiri_cluster_levels <- sort(unique(lparakefiri_cl))
lparakefiri_cluster_cols <- paletteer_d('ggthemes::Classic_10', n = length(lparakefiri_cluster_levels)) |>
  as.character()
names(lparakefiri_cluster_cols) <- as.character(lparakefiri_cluster_levels)

lparakefiri_tip_df <- tibble(
  label   = lparakefiri_tr$tip.label,
  cluster = factor(lparakefiri_cl[lparakefiri_tr$tip.label], levels = lparakefiri_cluster_levels)
)

# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

lparakefiri_site_cols <- c('Kabardino-Balkaria' = '#4FAFA8', 'Karachay-Cherkessia' = '#E8403D')

# Build TREE coordinates using ape (base), then extract from last_plot.phylo
# This gives us x/y positions without using ggtree/stat_tree.

plot(lparakefiri_tr, show.tip.label = FALSE)
lparakefiri_pp <- get('last_plot.phylo', envir = .PlotPhyloEnv)

# edge matrix: each row is (parent, child)
lparakefiri_edges <- as.data.frame(lparakefiri_tr$edge)
colnames(lparakefiri_edges) <- c('parent', 'child')

# Parent/child coordinates
lparakefiri_edges <- lparakefiri_edges |>
  mutate(
    x_parent = lparakefiri_pp$xx[parent],
    y_parent = lparakefiri_pp$yy[parent],
    x_child  = lparakefiri_pp$xx[child],
    y_child  = lparakefiri_pp$yy[child]
  )

# Build rectangular segments
lparakefiri_edges_h <- lparakefiri_edges |>
  transmute(x = x_parent, y = y_parent, xend = x_child, yend = y_parent)

lparakefiri_edges_v <- lparakefiri_edges |>
  transmute(x = x_child, y = y_parent, xend = x_child, yend = y_child)

# Tip coordinates
lparakefiri_tip_coords <- tibble(
  label = lparakefiri_tr$tip.label,
  node  = seq_along(lparakefiri_tr$tip.label),
  x     = lparakefiri_pp$xx[seq_along(lparakefiri_tr$tip.label)],
  y     = lparakefiri_pp$yy[seq_along(lparakefiri_tr$tip.label)]
) |>
  left_join(lparakefiri_tip_df, by = 'label') |>
  mutate(
    collection_site = if_else(label %in% kab_samples,
                              'Kabardino-Balkaria',
                              'Karachay-Cherkessia'),
    collection_site = factor(collection_site,
                             levels = c('Kabardino-Balkaria',
                                        'Karachay-Cherkessia'))
  )

# Order tips top-to-bottom
lparakefiri_tip_order <- lparakefiri_tip_coords |>
  arrange(desc(y)) |>
  mutate(tip_index = row_number()) |>
  select(label, tip_index)


# Plot the data -----------------------------------------------------------

lparakefiri_tip_label_offset <- 0.025 *
  diff(range(lparakefiri_edges$x_parent,
             lparakefiri_edges$x_child, na.rm = TRUE))

lparakefiri_tree_x_range <- diff(range(lparakefiri_edges$x_parent,
                                       lparakefiri_edges$x_child, na.rm = TRUE))

lparakefiri_leader_x_end <- max(lparakefiri_edges$x_child, na.rm = TRUE) +
  (0.12 * lparakefiri_tree_x_range)

lparakefiri_label_x <- lparakefiri_leader_x_end +
  (0.02 * lparakefiri_tree_x_range)

lparakefiri_dot_start_offset <- 0.04 * lparakefiri_tree_x_range

lparakefiri_tree_left_x <- min(lparakefiri_edges$x_parent,
                          lparakefiri_edges$x_child, na.rm = TRUE) -
  (0.05 * lparakefiri_tree_x_range)

# Tree segments
lparakefiri_tree_plot <- ggplot() +
  geom_segment(data = lparakefiri_edges_h,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.5) +
  geom_segment(data = lparakefiri_edges_v,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.5) +
  
  # Caucasus Region annotation
  geom_point(
    data = lparakefiri_tip_coords |> left_join(lparakefiri_tip_order, by = 'label'),
    aes(x = lparakefiri_tree_left_x, y = y,
        shape = collection_site, color = collection_site),
    size = 3,
    stroke = 0
  ) +
  scale_shape_manual(
    name   = 'Caucasus Region',
    values = c('Kabardino-Balkaria' = 17,
               'Karachay-Cherkessia'  = 15)
  ) +
  scale_color_manual(
    name   = 'Caucasus Region',
    values = lparakefiri_site_cols
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 4)),
    color = guide_legend(order = 1)
  ) +
  
  ggnewscale::new_scale_color() +
  
  # Cluster legend
  geom_point(
    data = lparakefiri_tip_coords |> left_join(lparakefiri_tip_order, by = 'label'),
    aes(x = x + lparakefiri_tip_label_offset * 0.3,
        y = y, color = cluster),
    size = 2
  ) +
  geom_segment(
    data = lparakefiri_tip_coords |> left_join(lparakefiri_tip_order, by = 'label'),
    aes(x = x + lparakefiri_dot_start_offset,
        y = y, xend = lparakefiri_leader_x_end,
        yend = y),
    linetype = 'dotted',
    linewidth = 0.4,
    color = 'grey70'
  ) +
  geom_text(
    data = lparakefiri_tip_coords |> left_join(lparakefiri_tip_order, by = 'label'),
    aes(x = lparakefiri_label_x, y = y,
        label = label, color = cluster),
    hjust = 0, size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = 'Clusters',
    values = lparakefiri_cluster_cols
  ) +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
    plot.margin = margin(0, 30, 0, 0),
    legend.title = element_text(face = 'bold', size = 12),
    legend.text  = element_text(size = 12)
  ) +
  guides(color = guide_legend(order = 2,
                              override.aes = list(shape = 16, size = 3)))

# Check the plot
lparakefiri_tree_plot


# Combine Tree and Phylogenetic Distance ----------------------------------

lparakefiri_combined_tree_dist_plot <-
  lparakefiri_tree_plot + lparakefiri_hm+
  plot_layout(widths = c(1, 1.2), guides = 'collect') +
  plot_annotation(
    title = '(B)',
    subtitle = expression(bolditalic('Lentilactobacillus parakefiri')),
    theme = theme(
      plot.title = element_text(face = 'bold', size = 12),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(face = 'bold', size = 12),
      legend.text  = element_text(size = 12)
    )
  )

# Check the plot
lparakefiri_combined_tree_dist_plot

# Save the plot
ggsave(
  lparakefiri_combined_tree_dist_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/lparakefiri_combined_plot.tiff',
  height = 15,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)