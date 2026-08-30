# StrainPhlAn t__SGB7221 Lentilactobacillus kefiri Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/tree_and_phylogenetic_distance
# 14/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkaria' = '#4FAFA8', 'Karachay-Cherkessia' = '#E8403D')


# Load the data -----------------------------------------------------------

lkefiri_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7221_Lentilactobacillus_kefiri_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
lkefiri_mat_df <- lkefiri_data |>
  bind_rows(
    lkefiri_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
lkefiri_all_ids <- sort(unique(c(lkefiri_mat_df$strain_id_1, lkefiri_mat_df$strain_id_2)))
lkefiri_diag_df <- tibble(
  strain_id_1 = lkefiri_all_ids,
  strain_id_2 = lkefiri_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
lkefiri_mat_df <- bind_rows(lkefiri_diag_df, lkefiri_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)

# Create region annotation table (one row per strain)
lkefiri_region_df <- tibble(strain_id = lkefiri_all_ids) |>
  mutate(
    collection_site = if_else(strain_id %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia'),
    collection_site = factor(collection_site, levels = c('Kabardino-Balkaria', 'Karachay-Cherkessia'))
  )

# Keep only sites that actually exist in this dataset (prevents KC showing in legend when absent)
lkefiri_sites_present <- lkefiri_region_df |>
  distinct(collection_site) |>
  pull(collection_site) |>
  as.character()

# Convert heatmap coordinates to numeric indices (so we can place annotations cleanly)
lkefiri_index_df <- tibble(
  strain_id = lkefiri_all_ids,
  x = seq_along(lkefiri_all_ids)
)

lkefiri_index_df_y <- tibble(
  strain_id = rev(lkefiri_all_ids),
  y = seq_along(rev(lkefiri_all_ids))
)

lkefiri_plot_df <- lkefiri_mat_df |>
  left_join(lkefiri_index_df, by = c('strain_id_1' = 'strain_id')) |>
  left_join(lkefiri_index_df_y, by = c('strain_id_2' = 'strain_id')) |>
  rename(x = x, y = y)

# Annotation positions (outside tile area, but NOT stretching axes)
lkefiri_n <- length(lkefiri_all_ids)

lkefiri_hm_left_x   <- 0.25  # left of first tile center (tile centers start at 1)
lkefiri_hm_bottom_y <- 0.25  # below the last tile row


# Plot the data -----------------------------------------------------------

lkefiri_hm <- ggplot(lkefiri_plot_df, aes(x = x, y = y, fill = phylogenetic_distance)) +
  geom_tile() +
  
  # Add distance values inside each tile
  geom_text(
    aes(label = sprintf('%.2f', phylogenetic_distance)),
    color = 'white',
    size  = 3.5
  ) +
  
  # LEFT annotation (Caucasus Region) aligned with rows
  geom_point(
    data = lkefiri_region_df |>
      left_join(lkefiri_index_df_y, by = 'strain_id'),
    aes(
      x = lkefiri_hm_left_x,
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
    data = lkefiri_region_df |>
      left_join(lkefiri_index_df, by = 'strain_id'),
    aes(
      x = x,
      y = lkefiri_hm_bottom_y,
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
    breaks = lkefiri_sites_present,
    drop   = TRUE
  ) +
  scale_color_manual(
    name   = 'Caucasus\nRegion',
    values = site_cols,
    breaks = lkefiri_sites_present,
    drop   = TRUE
  ) +
  
  # Heatmap fill scale
  scale_fill_gradient(
    name = 'Phylogenetic\nDistance',
    low  = '#081d3a',
    high = '#56b1f7',
    limits = c(0, max(lkefiri_plot_df$phylogenetic_distance, na.rm = TRUE)),
    na.value = 'white'
  ) +
  
  # Axis labels using numeric breaks -> strain IDs
  scale_x_continuous(
    breaks = seq_len(lkefiri_n),
    labels = lkefiri_all_ids,
    expand = c(0, 0),
    limits = c(0.0, lkefiri_n + 0.5)
  ) +
  scale_y_continuous(
    breaks = seq_len(lkefiri_n),
    labels = rev(lkefiri_all_ids),
    expand = c(0, 0),
    limits = c(lkefiri_n + 0.5, 0.0)
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
lkefiri_hm


# Making the Script Unique for L. kefiri ----------------------------------

# Load the necessary packages.

library(ape)
library(Biostrings)
library(paletteer)
library(patchwork)
library(ggnewscale)


# Load the data -----------------------------------------------------------

lkefiri_tree_file <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7221_Lentilactobacillus_kefiri/RAxML_bestTree.t__SGB7221.StrainPhlAn4.tre'
lkefiri_msa_file  <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7221_Lentilactobacillus_kefiri/t__SGB7221.StrainPhlAn4_concatenated.aln'

# Read tree + MSA
lkefiri_tr  <- read.tree(lkefiri_tree_file)
lkefiri_msa <- readDNAStringSet(lkefiri_msa_file, format = 'fasta') # width = 42

# Match/reorder
lkefiri_keep <- intersect(lkefiri_tr$tip.label, names(lkefiri_msa))
if (length(lkefiri_keep) == 0) stop('No overlap between tree tip labels and MSA names.')

lkefiri_tr  <- keep.tip(lkefiri_tr, lkefiri_keep)
lkefiri_msa <- lkefiri_msa[lkefiri_tr$tip.label]


# Define clusters (k <= 10 for Classic_10) --------------------------------

lkefiri_k <- 2
if (lkefiri_k > 10) stop('k must be <= 10 when using ggthemes::Classic_10.')

lkefiri_d  <- cophenetic.phylo(lkefiri_tr)
lkefiri_hc <- hclust(as.dist(lkefiri_d), method = 'average')
lkefiri_cl <- cutree(lkefiri_hc, k = lkefiri_k)
lkefiri_cl <- setNames(as.integer(lkefiri_cl), names(lkefiri_cl))

lkefiri_cluster_levels <- sort(unique(lkefiri_cl))
lkefiri_cluster_cols <- paletteer_d('ggthemes::Classic_10', n = length(lkefiri_cluster_levels)) |>
  as.character()
names(lkefiri_cluster_cols) <- as.character(lkefiri_cluster_levels)

lkefiri_tip_df <- tibble(
  label   = lkefiri_tr$tip.label,
  cluster = factor(lkefiri_cl[lkefiri_tr$tip.label], levels = lkefiri_cluster_levels)
)

# Build TREE coordinates using ape (base), then extract from last_plot.phylo
# This gives us x/y positions without using ggtree/stat_tree.

plot(lkefiri_tr, show.tip.label = FALSE)  # draws to device but we only need coordinates
lkefiri_pp <- get('last_plot.phylo', envir = .PlotPhyloEnv)

# edge matrix: each row is (parent, child)
lkefiri_edges <- as.data.frame(lkefiri_tr$edge)
colnames(lkefiri_edges) <- c('parent', 'child')

# Parent/child coordinates
lkefiri_edges <- lkefiri_edges |>
  mutate(
    x_parent = lkefiri_pp$xx[parent],
    y_parent = lkefiri_pp$yy[parent],
    x_child  = lkefiri_pp$xx[child],
    y_child  = lkefiri_pp$yy[child]
  )

# Build rectangular segments (horizontal + vertical) so the tree is standard rectangular
lkefiri_edges_h <- lkefiri_edges |>
  transmute(
    x    = x_parent,
    y    = y_parent,
    xend = x_child,
    yend = y_parent
  )

lkefiri_edges_v <- lkefiri_edges |>
  transmute(
    x    = x_child,
    y    = y_parent,
    xend = x_child,
    yend = y_child
  )

# Tip coordinates
lkefiri_tip_coords <- tibble(
  label = lkefiri_tr$tip.label,
  node  = seq_along(lkefiri_tr$tip.label),
  x     = lkefiri_pp$xx[seq_along(lkefiri_tr$tip.label)],
  y     = lkefiri_pp$yy[seq_along(lkefiri_tr$tip.label)]
) |>
  left_join(lkefiri_tip_df, by = 'label') |>
  mutate(
    collection_site = if_else(label %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia'),
    collection_site = factor(collection_site, levels = c('Kabardino-Balkaria', 'Karachay-Cherkessia'))
  )

# Order tips top-to-bottom as they appear in the plot (so MSA matches)
lkefiri_tip_order <- lkefiri_tip_coords |>
  arrange(desc(y)) |>
  mutate(tip_index = row_number()) |>
  select(label, tip_index)


# Plot the data -----------------------------------------------------------

# Tip label offset (move labels a bit farther from the tree line)
lkefiri_tip_label_offset <- 0.025 * diff(range(lkefiri_edges$x_parent, lkefiri_edges$x_child, na.rm = TRUE))

# Align dotted leaders and labels to fixed x-positions (adds a clean gap before the text)
lkefiri_tree_x_range <- diff(range(lkefiri_edges$x_parent, lkefiri_edges$x_child, na.rm = TRUE))
lkefiri_leader_x_end <- max(lkefiri_edges$x_child, na.rm = TRUE) + (0.12 * lkefiri_tree_x_range)  # where dots stop
lkefiri_label_x      <- lkefiri_leader_x_end + (0.02 * lkefiri_tree_x_range)                      # where text starts (gap)

# Offset so dotted line starts after the colored circle (in x units)
lkefiri_dot_start_offset <- 0.04 * lkefiri_tree_x_range   # adjust if needed

# Leader line length to align labels visually towards the MSA panel
lkefiri_leader_length <- 0.25 * diff(range(lkefiri_edges$x_parent, lkefiri_edges$x_child, na.rm = TRUE))

# Left annotation x position (Caucasus Region shapes)
lkefiri_tree_left_x <- min(lkefiri_edges$x_parent, lkefiri_edges$x_child, na.rm = TRUE) - (0.05 * lkefiri_tree_x_range)

# Tree segments
lkefiri_tree_plot <- ggplot() +
  geom_segment(
    data = lkefiri_edges_h,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  geom_segment(
    data = lkefiri_edges_v,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  
  # Caucasus Region annotation as shapes aligned to tips (triangle vs square)
  geom_point(
    data = lkefiri_tip_coords |> left_join(lkefiri_tip_order, by = 'label'),
    aes(x = lkefiri_tree_left_x, y = y, shape = collection_site, color = collection_site),
    size = 3,
    stroke = 0
  ) +
  scale_shape_manual(
    name   = 'Caucasus Region',
    values = c('Kabardino-Balkaria' = 17, 'Karachay-Cherkessia' = 15)
  ) +
  scale_color_manual(
    name   = 'Caucasus Region',
    values = site_cols
  ) +
  guides(
    shape = guide_legend(order = 1, override.aes = list(size = 4)),
    color = guide_legend(order = 1)
  ) +
  
  # start fresh color scale so clusters can have their own legend/colors
  ggnewscale::new_scale_color() +
  
  # Cluster legend as colored shapes (not letters)
  geom_point(
    data = lkefiri_tip_coords |> left_join(lkefiri_tip_order, by = 'label'),
    aes(x = x + lkefiri_tip_label_offset * 0.3, y = y, color = cluster),
    size = 2
  ) +
  # Dotted leader lines to visually connect tips to where the MSA starts
  geom_segment(
    data = lkefiri_tip_coords |> left_join(lkefiri_tip_order, by = 'label'),
    aes(
      x    = x + lkefiri_dot_start_offset,   # <-- start AFTER the circle
      y    = y,
      xend = lkefiri_leader_x_end,
      yend = y
    ),
    linetype = 'dotted',
    linewidth = 0.4,
    color = 'grey70'
  ) +
  geom_text(
    data = lkefiri_tip_coords |> left_join(lkefiri_tip_order, by = 'label'),
    aes(x = lkefiri_label_x, y = y, label = label, color = cluster),
    hjust = 0, size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = 'Clusters',
    values = lkefiri_cluster_cols,
    breaks = as.character(lkefiri_cluster_levels),
    labels = as.character(lkefiri_cluster_levels)
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
    legend.text = element_text(size = 12)
  ) +
  guides(
    # Put Clusters legend below Caucasus Region
    color = guide_legend(
      override.aes = list(shape = 16, size = 3),
      order = 2
    )
  )

# Check the plot
lkefiri_tree_plot


# Combine Tree and Phylogenetic Distance ----------------------------------

lkefiri_combined_tree_dist_plot <- lkefiri_tree_plot + lkefiri_hm +
  plot_layout(widths = c(1, 1.2), guides = 'collect') +
  plot_annotation(
    title = '(D)',
    subtitle = expression(bolditalic('Lentilactobacillus kefiri')),
    theme = theme(
      plot.title = element_text(face = 'bold', size = 12),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(face = 'bold', size = 12),
      legend.text  = element_text(size = 12)
    )
  )

# Check the plot
lkefiri_combined_tree_dist_plot

# Save the plot
ggsave(
  lkefiri_combined_tree_dist_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/lkefiri_combined_plot.tiff',
  height = 15,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine the Plots of Four Dominant Bacteria -----------------------------

combined_strain_plot <-
  (wrap_elements(full = edurans_combined_tree_dist_plot) /
     wrap_elements(full = llactis_combined_tree_dist_plot) /
     wrap_elements(full = lkefiranofaciens_combined_tree_dist_plot) /
     wrap_elements(full = lkefiri_combined_tree_dist_plot))

# Check the plot
combined_strain_plot

# Save the plot
ggsave(
  combined_strain_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/combined_strain_plot.tiff',
  height = 55,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)