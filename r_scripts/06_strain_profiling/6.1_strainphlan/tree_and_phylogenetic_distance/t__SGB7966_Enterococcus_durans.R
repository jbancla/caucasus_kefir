# StrainPhlAn t__SGB7966 Enterococcus durans Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/tree_and_phylogenetic_distance
# 14/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkaria' = '#4FAFA8', 'Karachay-Cherkessia' = '#E8403D')


# Load the data -----------------------------------------------------------

edurans_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7966_Enterococcus_durans_phylogenetic_distance_new.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
edurans_mat_df <- edurans_data |>
  bind_rows(
    edurans_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
edurans_all_ids <- sort(unique(c(edurans_mat_df$strain_id_1, edurans_mat_df$strain_id_2)))
edurans_diag_df <- tibble(
  strain_id_1 = edurans_all_ids,
  strain_id_2 = edurans_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
edurans_mat_df <- bind_rows(edurans_diag_df, edurans_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)

# Create region annotation table (one row per strain)
edurans_region_df <- tibble(strain_id = edurans_all_ids) |>
  mutate(
    collection_site = if_else(strain_id %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia'),
    collection_site = factor(collection_site, levels = c('Kabardino-Balkaria', 'Karachay-Cherkessia'))
  )

# Convert heatmap coordinates to numeric indices (so we can place annotations cleanly)
edurans_index_df <- tibble(
  strain_id = edurans_all_ids,
  x = seq_along(edurans_all_ids)
)

edurans_index_df_y <- tibble(
  strain_id = rev(edurans_all_ids),
  y = seq_along(rev(edurans_all_ids))
)

edurans_plot_df <- edurans_mat_df |>
  left_join(edurans_index_df, by = c('strain_id_1' = 'strain_id')) |>
  left_join(edurans_index_df_y, by = c('strain_id_2' = 'strain_id')) |>
  rename(x = x, y = y)

# Annotation positions (outside tile area, but NOT stretching axes)
edurans_n <- length(edurans_all_ids)

edurans_hm_left_x   <- 0.25  # left of first tile center (tile centers start at 1)
edurans_hm_bottom_y <- 0.25  # below the last tile row


# Plot the data -----------------------------------------------------------

edurans_hm <- ggplot(edurans_plot_df, aes(x = x, y = y, fill = phylogenetic_distance)) +
  geom_tile() +
  
  # Add distance values inside each tile
  geom_text(
    aes(label = sprintf('%.2f', phylogenetic_distance)),
    color = 'white',
    size  = 2.5 # size of phylogenetic distance value
  ) +
  
  # LEFT annotation (Caucasus Region) aligned with rows
  geom_point(
    data = edurans_region_df |>
      left_join(edurans_index_df_y |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = edurans_hm_left_x,
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
    data = edurans_region_df |>
      left_join(edurans_index_df |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = x,
      y = edurans_hm_bottom_y,
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
    values = c('Kabardino-Balkaria' = 17, 'Karachay-Cherkessia' = 15)
  ) +
  scale_color_manual(
    name   = 'Caucasus\nRegion',
    values = site_cols
  ) +
  
  # Heatmap fill scale
  scale_fill_gradient(
    name = 'Phylogenetic\nDistance',
    low  = '#081d3a',
    high = '#56b1f7',
    limits = c(0, max(edurans_plot_df$phylogenetic_distance, na.rm = TRUE)),
    na.value = 'white'
  ) +
  
  # Axis labels using numeric breaks -> strain IDs
  scale_x_continuous(
    breaks = seq_len(edurans_n),
    labels = edurans_all_ids,
    expand = c(0, 0),
    limits = c(0, edurans_n + 0.5)
  ) +
  scale_y_continuous(
    breaks = seq_len(edurans_n),
    labels = rev(edurans_all_ids),
    expand = c(0, 0),
    limits = c(edurans_n + 0.5, 0.0)
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
  guides(
    shape = 'none',
    color = 'none',
    fill  = guide_colorbar(order = 1)
  )

# Check the plot
edurans_hm


# Making the Script Unique for E. durans ----------------------------------

# Load the necessary packages.

library(ape)
library(Biostrings)
library(paletteer)
library(patchwork)
library(ggnewscale)


# Load the data -----------------------------------------------------------

edurans_tree_file <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7966_Enterococcus_durans/RAxML_bestTree.t__SGB7966.StrainPhlAn4.tre'
edurans_msa_file  <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7966_Enterococcus_durans/t__SGB7966.StrainPhlAn4_concatenated.aln'

# Read tree + MSA
edurans_tr  <- read.tree(edurans_tree_file)
edurans_msa <- readDNAStringSet(edurans_msa_file, format = 'fasta') # width = 213

# Match/reorder
edurans_keep <- intersect(edurans_tr$tip.label, names(edurans_msa))
if (length(edurans_keep) == 0) stop('No overlap between tree tip labels and MSA names.')

edurans_tr  <- keep.tip(edurans_tr, edurans_keep)
edurans_msa <- edurans_msa[edurans_tr$tip.label]


# Define clusters (k <= 10 for Classic_10) --------------------------------

edurans_k <- 3
if (edurans_k > 10) stop('k must be <= 10 when using ggthemes::Classic_10.')

edurans_d  <- cophenetic.phylo(edurans_tr)
edurans_hc <- hclust(as.dist(edurans_d), method = 'average')
edurans_cl <- cutree(edurans_hc, k = edurans_k)
edurans_cl <- setNames(as.integer(edurans_cl), names(edurans_cl))

edurans_cluster_levels <- sort(unique(edurans_cl))
edurans_cluster_cols <- paletteer_d('ggthemes::Classic_10', n = length(edurans_cluster_levels)) |>
  as.character()
names(edurans_cluster_cols) <- as.character(edurans_cluster_levels)

edurans_tip_df <- tibble(
  label   = edurans_tr$tip.label,
  cluster = factor(edurans_cl[edurans_tr$tip.label], levels = edurans_cluster_levels)
)

# Build TREE coordinates using ape (base), then extract from last_plot.phylo
# This gives us x/y positions without using ggtree/stat_tree.

plot(edurans_tr, show.tip.label = FALSE)  # draws to device but we only need coordinates
edurans_pp <- get('last_plot.phylo', envir = .PlotPhyloEnv)

# edge matrix: each row is (parent, child)
edurans_edges <- as.data.frame(edurans_tr$edge)
colnames(edurans_edges) <- c('parent', 'child')

# Parent/child coordinates
edurans_edges <- edurans_edges |>
  mutate(
    x_parent = edurans_pp$xx[parent],
    y_parent = edurans_pp$yy[parent],
    x_child  = edurans_pp$xx[child],
    y_child  = edurans_pp$yy[child]
  )

# Build rectangular segments (horizontal + vertical) so the tree is standard rectangular
edurans_edges_h <- edurans_edges |>
  transmute(
    x    = x_parent,
    y    = y_parent,
    xend = x_child,
    yend = y_parent
  )

edurans_edges_v <- edurans_edges |>
  transmute(
    x    = x_child,
    y    = y_parent,
    xend = x_child,
    yend = y_child
  )

# Tip coordinates
edurans_tip_coords <- tibble(
  label = edurans_tr$tip.label,
  node  = seq_along(edurans_tr$tip.label),
  x     = edurans_pp$xx[seq_along(edurans_tr$tip.label)],
  y     = edurans_pp$yy[seq_along(edurans_tr$tip.label)]
) |>
  left_join(edurans_tip_df, by = 'label') |>
  mutate(
    collection_site = if_else(label %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia'),
    collection_site = factor(collection_site, levels = c('Kabardino-Balkaria', 'Karachay-Cherkessia'))
  )

# Order tips top-to-bottom as they appear in the plot (so MSA matches)
edurans_tip_order <- edurans_tip_coords |>
  arrange(desc(y)) |>
  mutate(tip_index = row_number()) |>
  select(label, tip_index)


# Plot the data -----------------------------------------------------------

# Tip label offset (move labels a bit farther from the tree line)
edurans_tip_label_offset <- 0.025 * diff(range(edurans_edges$x_parent, edurans_edges$x_child, na.rm = TRUE))

# Align dotted leaders and labels to fixed x-positions (adds a clean gap before the text)
edurans_tree_x_range <- diff(range(edurans_edges$x_parent, edurans_edges$x_child, na.rm = TRUE))
edurans_leader_x_end <- max(edurans_edges$x_child, na.rm = TRUE) + (0.12 * edurans_tree_x_range)  # where dots stop
edurans_label_x      <- edurans_leader_x_end + (0.02 * edurans_tree_x_range)                      # where text starts (gap)

# Offset so dotted line starts after the colored circle (in x units)
edurans_dot_start_offset <- 0.04 * edurans_tree_x_range   # adjust if needed

# Leader line length to align labels visually towards the MSA panel
edurans_leader_length <- 0.25 * diff(range(edurans_edges$x_parent, edurans_edges$x_child, na.rm = TRUE))

# Left annotation x position (Caucasus Region shapes)
edurans_tree_left_x <- min(edurans_edges$x_parent, edurans_edges$x_child, na.rm = TRUE) - (0.05 * edurans_tree_x_range)

# Tree segments
edurans_tree_plot <- ggplot() +
  geom_segment(
    data = edurans_edges_h,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  geom_segment(
    data = edurans_edges_v,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  
  # Caucasus Region annotation as shapes aligned to tips (triangle vs square)
  geom_point(
    data = edurans_tip_coords |> left_join(edurans_tip_order, by = 'label'),
    aes(x = edurans_tree_left_x, y = y, shape = collection_site, color = collection_site),
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
    data = edurans_tip_coords |> left_join(edurans_tip_order, by = 'label'),
    aes(x = x + edurans_tip_label_offset * 0.3, y = y, color = cluster),
    size = 2
  ) +
  # Dotted leader lines to visually connect tips to where the MSA starts
  geom_segment(
    data = edurans_tip_coords |> left_join(edurans_tip_order, by = 'label'),
    aes(
      x    = x + edurans_dot_start_offset,   # <-- start AFTER the circle
      y    = y,
      xend = edurans_leader_x_end,
      yend = y
    ),
    linetype = 'dotted',
    linewidth = 0.4,
    color = 'grey70'
  ) +
  geom_text(
    data = edurans_tip_coords |> left_join(edurans_tip_order, by = 'label'),
    aes(x = edurans_label_x, y = y, label = label, color = cluster),
    hjust = 0, size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = 'Clusters',
    values = edurans_cluster_cols,
    breaks = as.character(edurans_cluster_levels),
    labels = as.character(edurans_cluster_levels)
  ) +
  coord_cartesian(clip = 'off') +
  theme_classic() +
  theme(
    plot.margin = margin(0, 30, 0, 0), # for spacing between tree and msa (activate to make white space between them)
    axis.line = element_blank(),
    axis.text = element_blank(),
    axis.ticks = element_blank(),
    axis.title = element_blank(),
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
edurans_tree_plot


# Combine Tree and Phylogenetic Distance ----------------------------------

edurans_combined_tree_dist_plot <- edurans_tree_plot + edurans_hm +
  plot_layout(widths = c(1, 1.2), guides = 'collect') +
  plot_annotation(
    title = '(A)',
    subtitle = expression(bolditalic('Enterococcus durans')),
    theme = theme(
      plot.title = element_text(face = 'bold', size = 12),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(face = 'bold', size = 12),
      legend.text  = element_text(size = 12)
    )
  )

# Check the plot
edurans_combined_tree_dist_plot

# Save the plot
ggsave(
  edurans_combined_tree_dist_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/edurans_combined_plot.tiff',
  height = 15,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)