# StrainPhlAn t__SGB7985 Lactococcus lactis Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/tree_and_phylogenetic_distance
# 14/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

site_cols <- c('Kabardino-Balkaria' = '#4FAFA8', 'Karachay-Cherkessia' = '#E8403D')


# Load the data -----------------------------------------------------------

llactis_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7985_Lactococcus_lactis_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
llactis_mat_df <- llactis_data |>
  bind_rows(
    llactis_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
llactis_all_ids <- sort(unique(c(llactis_mat_df$strain_id_1, llactis_mat_df$strain_id_2)))
llactis_diag_df <- tibble(
  strain_id_1 = llactis_all_ids,
  strain_id_2 = llactis_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
llactis_mat_df <- bind_rows(llactis_diag_df, llactis_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)

# Create region annotation table (one row per strain)
llactis_region_df <- tibble(strain_id = llactis_all_ids) |>
  mutate(
    collection_site = if_else(strain_id %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia'),
    collection_site = factor(collection_site, levels = c('Kabardino-Balkaria', 'Karachay-Cherkessia'))
  )

# Convert heatmap coordinates to numeric indices (so we can place annotations cleanly)
llactis_index_df <- tibble(
  strain_id = llactis_all_ids,
  x = seq_along(llactis_all_ids)
)

llactis_index_df_y <- tibble(
  strain_id = rev(llactis_all_ids),
  y = seq_along(rev(llactis_all_ids))
)

llactis_plot_df <- llactis_mat_df |>
  left_join(llactis_index_df, by = c('strain_id_1' = 'strain_id')) |>
  left_join(llactis_index_df_y, by = c('strain_id_2' = 'strain_id')) |>
  rename(x = x, y = y)

# Annotation positions (outside tile area, but NOT stretching axes)
llactis_n <- length(llactis_all_ids)

llactis_hm_left_x   <- 0.25  # left of first tile center (tile centers start at 1)
llactis_hm_bottom_y <- 0.25  # below the last tile row


# Plot the data -----------------------------------------------------------

llactis_hm <- ggplot(llactis_plot_df, aes(x = x, y = y, fill = phylogenetic_distance)) +
  geom_tile() +
  
  # Add distance values inside each tile
  geom_text(
    aes(label = sprintf('%.2f', phylogenetic_distance)),
    color = 'white',
    size  = 3
  ) +
  
  # LEFT annotation (Caucasus Region) aligned with rows
  geom_point(
    data = llactis_region_df |>
      left_join(llactis_index_df_y |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = llactis_hm_left_x,
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
    data = llactis_region_df |>
      left_join(llactis_index_df |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = x,
      y = llactis_hm_bottom_y,
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
    limits = c(0, max(llactis_plot_df$phylogenetic_distance, na.rm = TRUE)),
    na.value = 'white'
  ) +
  
  # Axis labels using numeric breaks -> strain IDs
  scale_x_continuous(
    breaks = seq_len(llactis_n),
    labels = llactis_all_ids,
    expand = c(0, 0),
    limits = c(0.0, llactis_n + 0.5)
  ) +
  scale_y_continuous(
    breaks = seq_len(llactis_n),
    labels = rev(llactis_all_ids),
    expand = c(0, 0),
    limits = c(llactis_n + 0.5, 0.0)
  ) +
  coord_equal() +
  theme_classic() +
  theme(
    axis.title   = element_blank(),
    axis.text   = element_text(size = 12),
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
llactis_hm


# Making the Script Unique for L. lactis ----------------------------------

# Load the necessary packages.

library(ape)
library(Biostrings)
library(paletteer)
library(patchwork)
library(ggnewscale)


# Load the data -----------------------------------------------------------

llactis_tree_file <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7985_Lactococcus_lactis/RAxML_bestTree.t__SGB7985.StrainPhlAn4.tre'
llactis_msa_file  <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7985_Lactococcus_lactis/t__SGB7985.StrainPhlAn4_concatenated.aln'

# Read tree + MSA
llactis_tr  <- read.tree(llactis_tree_file)
llactis_msa <- readDNAStringSet(llactis_msa_file, format = 'fasta') # width = 53

# Match/reorder
llactis_keep <- intersect(llactis_tr$tip.label, names(llactis_msa))
if (length(llactis_keep) == 0) stop('No overlap between tree tip labels and MSA names.')

llactis_tr  <- keep.tip(llactis_tr, llactis_keep)
llactis_msa <- llactis_msa[llactis_tr$tip.label]


# Define clusters (k <= 10 for Classic_10) --------------------------------

llactis_k <- 3
if (llactis_k > 10) stop('k must be <= 10 when using ggthemes::Classic_10.')

llactis_d  <- cophenetic.phylo(llactis_tr)
llactis_hc <- hclust(as.dist(llactis_d), method = 'average')
llactis_cl <- cutree(llactis_hc, k = llactis_k)
llactis_cl <- setNames(as.integer(llactis_cl), names(llactis_cl))

llactis_cluster_levels <- sort(unique(llactis_cl))
llactis_cluster_cols <- paletteer_d('ggthemes::Classic_10', n = length(llactis_cluster_levels)) |>
  as.character()
names(llactis_cluster_cols) <- as.character(llactis_cluster_levels)

llactis_tip_df <- tibble(
  label   = llactis_tr$tip.label,
  cluster = factor(llactis_cl[llactis_tr$tip.label], levels = llactis_cluster_levels)
)

# Build TREE coordinates using ape (base), then extract from last_plot.phylo
# This gives us x/y positions without using ggtree/stat_tree.

plot(llactis_tr, show.tip.label = FALSE)  # draws to device but we only need coordinates
llactis_pp <- get('last_plot.phylo', envir = .PlotPhyloEnv)

# edge matrix: each row is (parent, child)
llactis_edges <- as.data.frame(llactis_tr$edge)
colnames(llactis_edges) <- c('parent', 'child')

# Parent/child coordinates
llactis_edges <- llactis_edges |>
  mutate(
    x_parent = llactis_pp$xx[parent],
    y_parent = llactis_pp$yy[parent],
    x_child  = llactis_pp$xx[child],
    y_child  = llactis_pp$yy[child]
  )

# Build rectangular segments (horizontal + vertical) so the tree is standard rectangular
llactis_edges_h <- llactis_edges |>
  transmute(
    x    = x_parent,
    y    = y_parent,
    xend = x_child,
    yend = y_parent
  )

llactis_edges_v <- llactis_edges |>
  transmute(
    x    = x_child,
    y    = y_parent,
    xend = x_child,
    yend = y_child
  )

# Tip coordinates
llactis_tip_coords <- tibble(
  label = llactis_tr$tip.label,
  node  = seq_along(llactis_tr$tip.label),
  x     = llactis_pp$xx[seq_along(llactis_tr$tip.label)],
  y     = llactis_pp$yy[seq_along(llactis_tr$tip.label)]
) |>
  left_join(llactis_tip_df, by = 'label') |>
  mutate(
    collection_site = if_else(label %in% kab_samples, 'Kabardino-Balkaria', 'Karachay-Cherkessia'),
    collection_site = factor(collection_site, levels = c('Kabardino-Balkaria', 'Karachay-Cherkessia'))
  )

# Order tips top-to-bottom as they appear in the plot (so MSA matches)
llactis_tip_order <- llactis_tip_coords |>
  arrange(desc(y)) |>
  mutate(tip_index = row_number()) |>
  select(label, tip_index)

# Plot the data -----------------------------------------------------------

# Tip label offset (move labels a bit farther from the tree line)
llactis_tip_label_offset <- 0.025 * diff(range(llactis_edges$x_parent, llactis_edges$x_child, na.rm = TRUE))

# Align dotted leaders and labels to fixed x-positions (adds a clean gap before the text)
llactis_tree_x_range <- diff(range(llactis_edges$x_parent, llactis_edges$x_child, na.rm = TRUE))
llactis_leader_x_end <- max(llactis_edges$x_child, na.rm = TRUE) + (0.12 * llactis_tree_x_range)  # where dots stop
llactis_label_x      <- llactis_leader_x_end + (0.02 * llactis_tree_x_range)                      # where text starts (gap)

# Offset so dotted line starts after the colored circle (in x units)
llactis_dot_start_offset <- 0.04 * llactis_tree_x_range   # adjust if needed

# Leader line length to align labels visually towards the MSA panel
llactis_leader_length <- 0.25 * diff(range(llactis_edges$x_parent, llactis_edges$x_child, na.rm = TRUE))

# Left annotation x position (Caucasus Region shapes)
llactis_tree_left_x <- min(llactis_edges$x_parent, llactis_edges$x_child, na.rm = TRUE) - (0.05 * llactis_tree_x_range)

# Tree segments
llactis_tree_plot <- ggplot() +
  geom_segment(
    data = llactis_edges_h,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  geom_segment(
    data = llactis_edges_v,
    aes(x = x, y = y, xend = xend, yend = yend),
    linewidth = 0.5
  ) +
  
  # Caucasus Region annotation as shapes aligned to tips (triangle vs square)
  geom_point(
    data = llactis_tip_coords |> left_join(llactis_tip_order, by = 'label'),
    aes(x = llactis_tree_left_x, y = y, shape = collection_site, color = collection_site),
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
    data = llactis_tip_coords |> left_join(llactis_tip_order, by = 'label'),
    aes(x = x + llactis_tip_label_offset * 0.3, y = y, color = cluster),
    size = 2
  ) +
  # Dotted leader lines to visually connect tips to where the MSA starts
  geom_segment(
    data = llactis_tip_coords |> left_join(llactis_tip_order, by = 'label'),
    aes(
      x    = x + llactis_dot_start_offset,   # <-- start AFTER the circle
      y    = y,
      xend = llactis_leader_x_end,
      yend = y
    ),
    linetype = 'dotted',
    linewidth = 0.4,
    color = 'grey70'
  ) +
  geom_text(
    data = llactis_tip_coords |> left_join(llactis_tip_order, by = 'label'),
    aes(x = llactis_label_x, y = y, label = label, color = cluster),
    hjust = 0, size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = 'Clusters',
    values = llactis_cluster_cols,
    breaks = as.character(llactis_cluster_levels),
    labels = as.character(llactis_cluster_levels)
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
llactis_tree_plot


# Combine Tree and Phylogenetic Distance ----------------------------------

llactis_combined_tree_dist_plot <- llactis_tree_plot + llactis_hm +
  plot_layout(widths = c(1, 1.2), guides = 'collect') +
  plot_annotation(
    title = '(B)',
    subtitle = expression(bolditalic('Lactococcus lactis')),
    theme = theme(
      plot.title = element_text(face = 'bold', size = 12),
      plot.subtitle = element_text(size  = 12, hjust = 0.5),
      legend.title = element_text(face = 'bold', size = 12),
      legend.text  = element_text(size = 12)
    )
  )

# Check the plot
llactis_combined_tree_dist_plot

# Save the plot
ggsave(
  llactis_combined_tree_dist_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/llactis_combined_plot.tiff',
  height = 15,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)