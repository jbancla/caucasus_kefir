# StrainPhlAn t__SGB7046 Lactobacillus kefiranofaciens Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/tree_and_phylogenetic_distance
# 14/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkaria' = '#4FAFA8', 'Karachay-Cherkessia' = '#E8403D')


# Load the data -----------------------------------------------------------

lkefiranofaciens_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7046_Lactobacillus_kefiranofaciens_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
lkefiranofaciens_mat_df <- lkefiranofaciens_data |>
  bind_rows(
    lkefiranofaciens_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
lkefiranofaciens_all_ids <- sort(unique(c(lkefiranofaciens_mat_df$strain_id_1, lkefiranofaciens_mat_df$strain_id_2)))
lkefiranofaciens_diag_df <- tibble(
  strain_id_1 = lkefiranofaciens_all_ids,
  strain_id_2 = lkefiranofaciens_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
lkefiranofaciens_mat_df <- bind_rows(lkefiranofaciens_diag_df, lkefiranofaciens_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)

# Create region annotation table (one row per strain)
lkefiranofaciens_region_df <- tibble(strain_id = lkefiranofaciens_all_ids) |>
  mutate(
    collection_site = if_else(strain_id %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia'),
    collection_site = factor(collection_site, levels = c('Kabardino-Balkaria', 'Karachay-Cherkessia'))
  )

# Keep only sites that actually exist in this dataset (prevents KC showing in legend when absent)
lkefiranofaciens_sites_present <- lkefiranofaciens_region_df |>
  distinct(collection_site) |>
  pull(collection_site) |>
  as.character()

# Convert heatmap coordinates to numeric indices (so we can place annotations cleanly)
lkefiranofaciens_index_df <- tibble(
  strain_id = lkefiranofaciens_all_ids,
  x = seq_along(lkefiranofaciens_all_ids)
)

lkefiranofaciens_index_df_y <- tibble(
  strain_id = rev(lkefiranofaciens_all_ids),
  y = seq_along(rev(lkefiranofaciens_all_ids))
)

lkefiranofaciens_plot_df <- lkefiranofaciens_mat_df |>
  left_join(lkefiranofaciens_index_df, by = c('strain_id_1' = 'strain_id')) |>
  left_join(lkefiranofaciens_index_df_y, by = c('strain_id_2' = 'strain_id')) |>
  rename(x = x, y = y)

# Annotation positions (outside tile area, but NOT stretching axes)
lkefiranofaciens_n <- length(lkefiranofaciens_all_ids)

lkefiranofaciens_hm_left_x   <- 0.25  # left of first tile center (tile centers start at 1)
lkefiranofaciens_hm_bottom_y <- 0.25  # below the last tile row


# Plot the data -----------------------------------------------------------

lkefiranofaciens_hm <- ggplot(lkefiranofaciens_plot_df, aes(x = x, y = y, fill = phylogenetic_distance)) +
  geom_tile() +
  
  # Add distance values inside each tile
  geom_text(
    aes(label = sprintf('%.2f', phylogenetic_distance)),
    color = 'white',
    size  = 3.5
  ) +
  
  # LEFT annotation (Caucasus Region) aligned with rows
  geom_point(
    data = lkefiranofaciens_region_df |>
      left_join(lkefiranofaciens_index_df_y, by = 'strain_id'),
    aes(
      x = lkefiranofaciens_hm_left_x,
      y = y,
      shape = collection_site,
      color = collection_site
    ),
    inherit.aes = FALSE,
    size = 3.5,
    stroke = 0
  ) +
  
  # BOTTOM annotation (Caucasus Region) aligned with columns
  geom_point(
    data = lkefiranofaciens_region_df |>
      left_join(lkefiranofaciens_index_df, by = 'strain_id'),
    aes(
      x = x,
      y = lkefiranofaciens_hm_bottom_y,
      shape = collection_site,
      color = collection_site
    ),
    inherit.aes = FALSE,
    size = 3.5,
    stroke = 0
  ) +
  
  # Shapes and colors for Caucasus Region
  scale_shape_manual(
    name   = 'Caucasus\nRegion',
    values = c('Kabardino-Balkaria' = 17, 'Karachay-Cherkessia' = 15),
    breaks = lkefiranofaciens_sites_present,
    drop   = TRUE
  ) +
  scale_color_manual(
    name   = 'Caucasus\nRegion',
    values = site_cols,
    breaks = lkefiranofaciens_sites_present,
    drop   = TRUE
  ) +
  
  # Heatmap fill scale
  scale_fill_gradient(
    name = 'Phylogenetic\nDistance',
    low  = '#081d3a',
    high = '#56b1f7',
    limits = c(0, max(lkefiranofaciens_plot_df$phylogenetic_distance, na.rm = TRUE)),
    na.value = 'white'
  ) +
  
  # Axis labels using numeric breaks -> strain IDs
  scale_x_continuous(
    breaks = seq_len(lkefiranofaciens_n),
    labels = lkefiranofaciens_all_ids,
    expand = c(0, 0),
    limits = c(0.0, lkefiranofaciens_n + 0.5)
  ) +
  scale_y_continuous(
    breaks = seq_len(lkefiranofaciens_n),
    labels = rev(lkefiranofaciens_all_ids),
    expand = c(0, 0),
    limits = c(lkefiranofaciens_n + 0.5, 0.0)
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    axis.title   = element_blank(),
    axis.text = element_text(size = 12),
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
lkefiranofaciens_hm


# Making the Script Unique for L. kefiranofaciens -------------------------

# Load the necessary packages.

library(ape)
library(Biostrings)
library(paletteer)
library(patchwork)
library(ggnewscale)


# Load the data -----------------------------------------------------------

lkefiranofaciens_tree_file <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7046_Lactobacillus_kefiranofaciens/RAxML_bestTree.t__SGB7046.StrainPhlAn4.tre'
lkefiranofaciens_msa_file  <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7046_Lactobacillus_kefiranofaciens/t__SGB7046.StrainPhlAn4_concatenated.aln'

# Read tree + MSA
lkefiranofaciens_tr  <- read.tree(lkefiranofaciens_tree_file)
lkefiranofaciens_msa <- readDNAStringSet(lkefiranofaciens_msa_file, format = 'fasta') # width = 98

# Match/reorder
lkefiranofaciens_keep <- intersect(lkefiranofaciens_tr$tip.label, names(lkefiranofaciens_msa))
if (length(lkefiranofaciens_keep) == 0) stop('No overlap between tree tip labels and MSA names.')

lkefiranofaciens_tr  <- keep.tip(lkefiranofaciens_tr, lkefiranofaciens_keep)
lkefiranofaciens_msa <- lkefiranofaciens_msa[lkefiranofaciens_tr$tip.label]


# Define clusters (k <= 10 for Classic_10) --------------------------------

lkefiranofaciens_k <- 3
if (lkefiranofaciens_k > 10) stop('k must be <= 10 when using ggthemes::Classic_10.')

lkefiranofaciens_d  <- cophenetic.phylo(lkefiranofaciens_tr)
lkefiranofaciens_hc <- hclust(as.dist(lkefiranofaciens_d), method = 'average')
lkefiranofaciens_cl <- cutree(lkefiranofaciens_hc, k = lkefiranofaciens_k)
lkefiranofaciens_cl <- setNames(as.integer(lkefiranofaciens_cl), names(lkefiranofaciens_cl))

lkefiranofaciens_cluster_levels <- sort(unique(lkefiranofaciens_cl))
lkefiranofaciens_cluster_cols <- paletteer_d('ggthemes::Classic_10',
                                             n = length(lkefiranofaciens_cluster_levels)) |>
  as.character()
names(lkefiranofaciens_cluster_cols) <- as.character(lkefiranofaciens_cluster_levels)

lkefiranofaciens_tip_df <- tibble(
  label   = lkefiranofaciens_tr$tip.label,
  cluster = factor(lkefiranofaciens_cl[lkefiranofaciens_tr$tip.label],
                   levels = lkefiranofaciens_cluster_levels)
)

# Build TREE coordinates using ape (base), then extract from last_plot.phylo

plot(lkefiranofaciens_tr, show.tip.label = FALSE)
lkefiranofaciens_pp <- get('last_plot.phylo', envir = .PlotPhyloEnv)

lkefiranofaciens_edges <- as.data.frame(lkefiranofaciens_tr$edge)
colnames(lkefiranofaciens_edges) <- c('parent', 'child')

lkefiranofaciens_edges <- lkefiranofaciens_edges |>
  mutate(
    x_parent = lkefiranofaciens_pp$xx[parent],
    y_parent = lkefiranofaciens_pp$yy[parent],
    x_child  = lkefiranofaciens_pp$xx[child],
    y_child  = lkefiranofaciens_pp$yy[child]
  )

lkefiranofaciens_edges_h <- lkefiranofaciens_edges |>
  transmute(x = x_parent, y = y_parent, xend = x_child, yend = y_parent)

lkefiranofaciens_edges_v <- lkefiranofaciens_edges |>
  transmute(x = x_child, y = y_parent, xend = x_child, yend = y_child)

# Tip coordinates
lkefiranofaciens_tip_coords <- tibble(
  label = lkefiranofaciens_tr$tip.label,
  node  = seq_along(lkefiranofaciens_tr$tip.label),
  x     = lkefiranofaciens_pp$xx[seq_along(lkefiranofaciens_tr$tip.label)],
  y     = lkefiranofaciens_pp$yy[seq_along(lkefiranofaciens_tr$tip.label)]
) |>
  left_join(lkefiranofaciens_tip_df, by = 'label') |>
  mutate(
    collection_site = if_else(label %in% kab_samples,
                              'Kabardino-Balkaria',
                              'Karachay-Cherkessia'),
    collection_site = factor(collection_site,
                             levels = c('Kabardino-Balkaria',
                                        'Karachay-Cherkessia'))
  )

# Order tips top-to-bottom as they appear in the plot (so MSA matches)
lkefiranofaciens_tip_order <- lkefiranofaciens_tip_coords |>
  arrange(desc(y)) |>
  mutate(tip_index = row_number()) |>
  select(label, tip_index)


# Plot the data -----------------------------------------------------------

lkefiranofaciens_tip_label_offset <- 0.025 *
  diff(range(lkefiranofaciens_edges$x_parent,
             lkefiranofaciens_edges$x_child,
             na.rm = TRUE))

lkefiranofaciens_tree_x_range <- diff(range(lkefiranofaciens_edges$x_parent,
                                            lkefiranofaciens_edges$x_child,
                                            na.rm = TRUE))

lkefiranofaciens_leader_x_end <- max(lkefiranofaciens_edges$x_child,
                                     na.rm = TRUE) +
  (0.12 * lkefiranofaciens_tree_x_range)

lkefiranofaciens_label_x <- lkefiranofaciens_leader_x_end +
  (0.02 * lkefiranofaciens_tree_x_range)

lkefiranofaciens_dot_start_offset <- 0.04 *
  lkefiranofaciens_tree_x_range

lkefiranofaciens_leader_length <- 0.25 *
  diff(range(lkefiranofaciens_edges$x_parent,
             lkefiranofaciens_edges$x_child,
             na.rm = TRUE))

lkefiranofaciens_tree_left_x <- min(lkefiranofaciens_edges$x_parent,
                               lkefiranofaciens_edges$x_child,
                               na.rm = TRUE) -
  (0.05 * lkefiranofaciens_tree_x_range)


# Tree segments
lkefiranofaciens_tree_plot <- ggplot() +
  geom_segment(
    data = lkefiranofaciens_edges_h,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  geom_segment(
    data = lkefiranofaciens_edges_v,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  
  # Caucasus Region annotation
  geom_point(
    data = lkefiranofaciens_tip_coords,
    aes(x = lkefiranofaciens_tree_left_x,
        y = y,
        shape = collection_site,
        color = collection_site),
    size = 3
  ) +
  scale_shape_manual(
    name   = 'Caucasus Region',
    values = c('Kabardino-Balkaria' = 17,
               'Karachay-Cherkessia'    = 15)
  ) +
  scale_color_manual(
    name   = 'Caucasus Region',
    values = site_cols
  ) +
  guides(
    shape = guide_legend(order = 1),
    color = guide_legend(order = 1)
  ) +
  
  ggnewscale::new_scale_color() +
  
  # Cluster legend
  geom_point(
    data = lkefiranofaciens_tip_coords,
    aes(x = x + lkefiranofaciens_tip_label_offset * 0.3,
        y = y,
        color = cluster),
    size = 2
  ) +
  geom_segment(
    data = lkefiranofaciens_tip_coords,
    aes(
      x    = x + lkefiranofaciens_dot_start_offset,
      y    = y,
      xend = lkefiranofaciens_leader_x_end,
      yend = y
    ),
    linetype = 'dotted',
    linewidth = 0.4,
    color = 'grey70'
  ) +
  geom_text(
    data = lkefiranofaciens_tip_coords,
    aes(x = lkefiranofaciens_label_x,
        y = y,
        label = label,
        color = cluster),
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = 'Clusters',
    values = lkefiranofaciens_cluster_cols
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
  guides(
    color = guide_legend(
      override.aes = list(shape = 16, size = 3),
      order = 2
    )
  )

# Check the plot
lkefiranofaciens_tree_plot


# Combine Tree and Phylogenetic Distance ----------------------------------

lkefiranofaciens_combined_tree_dist_plot <-
  lkefiranofaciens_tree_plot + lkefiranofaciens_hm +
  plot_layout(widths = c(1, 1.2), guides = 'collect') +
  plot_annotation(
    title = '(C)',
    subtitle = expression(bolditalic('Lactobacillus kefiranofaciens')),
    theme = theme(
      plot.title = element_text(face = 'bold', size = 12),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(face = 'bold', size = 12),
      legend.text  = element_text(size = 12)
    )
  )

# Check the plot
lkefiranofaciens_combined_tree_dist_plot

# Save the plot
ggsave(
  lkefiranofaciens_combined_tree_dist_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/lkefiranofaciens_combined_plot.tiff',
  height = 15,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)
