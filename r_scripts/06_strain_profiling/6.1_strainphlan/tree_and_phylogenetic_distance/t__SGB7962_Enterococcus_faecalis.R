# StrainPhlAn t__SGB7962 Enterococcus faecalis Phylogenetic Distance
# ~/Documents/RStudio/mk_caucasus/scripts/06_strain_profiling/6.1_strainphlan/tree_and_phylogenetic_distance
# 17/02/2026

# Load the necessary packages.

library(tidyverse)


# Adding collection sites layer -------------------------------------------
# NOTE: E. faecalis has ONLY KC data here (so everything is Karachay-Cherkess)

site_cols <- c('Kabardino-Balkarian' = '#4FAFA8', 'Karachay-Cherkess' = '#E8403D')


# Load the data -----------------------------------------------------------

efaecalis_data <- read_tsv(
  'data/06_strain_profiling/6.1_strainphlan/6.1.5_phylogenetic_distances/t__SGB7962_Enterococcus_faecalis_phylogenetic_distance.tsv',
  show_col_types = FALSE
)

# Build symmetric matrix
efaecalis_mat_df <- efaecalis_data |>
  bind_rows(
    efaecalis_data |>
      transmute(
        strain_id_1 = strain_id_2,
        strain_id_2 = strain_id_1,
        phylogenetic_distance
      )
  ) |>
  distinct() |>
  mutate(phylogenetic_distance = as.numeric(phylogenetic_distance))

# Add diagonal = 0
efaecalis_all_ids <- sort(unique(c(efaecalis_mat_df$strain_id_1, efaecalis_mat_df$strain_id_2)))

efaecalis_diag_df <- tibble(
  strain_id_1 = efaecalis_all_ids,
  strain_id_2 = efaecalis_all_ids,
  phylogenetic_distance = 0
)

# IMPORTANT: bind diag_df FIRST so its 0-values win, then keep the first occurrence
efaecalis_mat_df <- bind_rows(efaecalis_diag_df, efaecalis_mat_df) |>
  distinct(strain_id_1, strain_id_2, .keep_all = TRUE)


# Heatmap -----------------------------------------------------------------

# Create region annotation table (one row per strain)
efaecalis_region_df <- tibble(strain_id = efaecalis_all_ids) |>
  mutate(
    collection_site = 'Karachay-Cherkess',
    collection_site = factor(
      collection_site,
      levels = c('Kabardino-Balkarian', 'Karachay-Cherkess')
    )
  )

# Keep only sites that actually exist in this dataset (prevents KB showing in legend when absent)
efaecalis_sites_present <- efaecalis_region_df |>
  distinct(collection_site) |>
  pull(collection_site) |>
  as.character()

# Convert heatmap coordinates to numeric indices (so we can place annotations cleanly)
efaecalis_index_df <- tibble(
  strain_id = efaecalis_all_ids,
  x = seq_along(efaecalis_all_ids)
)

efaecalis_index_df_y <- tibble(
  strain_id = rev(efaecalis_all_ids),
  y = seq_along(rev(efaecalis_all_ids))
)

efaecalis_plot_df <- efaecalis_mat_df |>
  left_join(efaecalis_index_df, by = c('strain_id_1' = 'strain_id')) |>
  left_join(efaecalis_index_df_y, by = c('strain_id_2' = 'strain_id')) |>
  rename(x = x, y = y)

# Annotation positions (outside tile area, but NOT stretching axes)
efaecalis_n <- length(efaecalis_all_ids)

efaecalis_hm_left_x   <- 0.25
efaecalis_hm_bottom_y <- 0.25


# Plot the data -----------------------------------------------------------

efaecalis_hm <- ggplot(efaecalis_plot_df, aes(x = x, y = y, fill = phylogenetic_distance)) +
  geom_tile() +
  
  # Add distance values inside each tile
  geom_text(
    aes(label = sprintf('%.2f', phylogenetic_distance)),
    color = 'white',
    size  = 3.5
  ) +
  
  # LEFT annotation (Caucasus Region) aligned with rows
  geom_point(
    data = efaecalis_region_df |>
      left_join(efaecalis_index_df_y |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = efaecalis_hm_left_x,
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
    data = efaecalis_region_df |>
      left_join(efaecalis_index_df |> rename(strain_id = strain_id), by = 'strain_id'),
    aes(
      x = x,
      y = efaecalis_hm_bottom_y,
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
    breaks = efaecalis_sites_present,
    drop   = TRUE
  ) +
  scale_color_manual(
    name   = 'Caucasus\nRegion',
    values = site_cols,
    breaks = efaecalis_sites_present,
    drop   = TRUE
  ) +
  
  # Heatmap fill scale
  scale_fill_gradient(
    name = 'Phylogenetic\nDistance',
    low  = '#081d3a',
    high = '#56b1f7',
    limits = c(0, max(efaecalis_plot_df$phylogenetic_distance, na.rm = TRUE)),
    na.value = 'white'
  ) +
  
  # Axis labels using numeric breaks -> strain IDs
  scale_x_continuous(
    breaks = seq_len(efaecalis_n),
    labels = efaecalis_all_ids,
    expand = c(0, 0),
    limits = c(0.0, efaecalis_n + 0.5)
  ) +
  scale_y_continuous(
    breaks = seq_len(efaecalis_n),
    labels = rev(efaecalis_all_ids),
    expand = c(0, 0),
    limits = c(efaecalis_n + 0.5, 0.0)
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
  
  guides(
    shape = 'none',
    color = 'none',
    fill  = guide_colorbar(order = 1)
  )

# Check the plot
efaecalis_hm


# Making the Script Unique for E. faecalis --------------------------------

# Load the necessary packages.

library(ape)
library(Biostrings)
library(paletteer)
library(patchwork)
library(ggnewscale)


# Load the data -----------------------------------------------------------

efaecalis_tree_file <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7962_Enterococcus_faecalis/RAxML_bestTree.t__SGB7962.StrainPhlAn4.tre'
efaecalis_msa_file  <- 'data/06_strain_profiling/6.1_strainphlan/6.1.4_output/t__SGB7962_Enterococcus_faecalis/t__SGB7962.StrainPhlAn4_concatenated.aln'

# Read tree + MSA
efaecalis_tr  <- read.tree(efaecalis_tree_file)
efaecalis_msa <- readDNAStringSet(efaecalis_msa_file, format = 'fasta') # width = 10

# Match/reorder
efaecalis_keep <- intersect(efaecalis_tr$tip.label, names(efaecalis_msa))
if (length(efaecalis_keep) == 0) stop('No overlap between tree tip labels and MSA names.')

efaecalis_tr  <- keep.tip(efaecalis_tr, efaecalis_keep)
efaecalis_msa <- efaecalis_msa[efaecalis_tr$tip.label]


# Define clusters (k <= 10 for Classic_10) --------------------------------

efaecalis_k <- 3
if (efaecalis_k > 10) stop('k must be <= 10 when using ggthemes::Classic_10.')

efaecalis_d  <- cophenetic.phylo(efaecalis_tr)
efaecalis_hc <- hclust(as.dist(efaecalis_d), method = 'average')
efaecalis_cl <- cutree(efaecalis_hc, k = efaecalis_k)
efaecalis_cl <- setNames(as.integer(efaecalis_cl), names(efaecalis_cl))

efaecalis_cluster_levels <- sort(unique(efaecalis_cl))
efaecalis_cluster_cols <- paletteer_d('ggthemes::Classic_10',
                                      n = length(efaecalis_cluster_levels)) |>
  as.character()
names(efaecalis_cluster_cols) <- as.character(efaecalis_cluster_levels)

efaecalis_tip_df <- tibble(
  label   = efaecalis_tr$tip.label,
  cluster = factor(efaecalis_cl[efaecalis_tr$tip.label],
                   levels = efaecalis_cluster_levels)
)


# Adding collection sites layer -------------------------------------------

kab_samples <- c('k1_m','k2_m','k3_m','k5_m','k6_m','k7_m','k8_m','k9_m',
                 'k1_g','k3_g','k5_g','k6_g','k7_g','k9_g')

efaecalis_site_cols <- c('Kabardino-Balkaria' = '#4FAFA8', 'Karachay-Cherkessia' = '#E8403D')

# Build TREE coordinates using ape (base), then extract from last_plot.phylo

plot(efaecalis_tr, show.tip.label = FALSE)
efaecalis_pp <- get('last_plot.phylo', envir = .PlotPhyloEnv)

# edge matrix
efaecalis_edges <- as.data.frame(efaecalis_tr$edge)
colnames(efaecalis_edges) <- c('parent', 'child')

# Parent/child coordinates
efaecalis_edges <- efaecalis_edges |>
  mutate(
    x_parent = efaecalis_pp$xx[parent],
    y_parent = efaecalis_pp$yy[parent],
    x_child  = efaecalis_pp$xx[child],
    y_child  = efaecalis_pp$yy[child]
  )

# Rectangular segments
efaecalis_edges_h <- efaecalis_edges |>
  transmute(x = x_parent, y = y_parent,
            xend = x_child, yend = y_parent)

efaecalis_edges_v <- efaecalis_edges |>
  transmute(x = x_child, y = y_parent,
            xend = x_child, yend = y_child)

# Tip coordinates
efaecalis_tip_coords <- tibble(
  label = efaecalis_tr$tip.label,
  node  = seq_along(efaecalis_tr$tip.label),
  x     = efaecalis_pp$xx[seq_along(efaecalis_tr$tip.label)],
  y     = efaecalis_pp$yy[seq_along(efaecalis_tr$tip.label)]
) |>
  left_join(efaecalis_tip_df, by = 'label') |>
  mutate(
    collection_site = if_else(label %in% kab_samples,
                              'Kabardino-Balkaria',
                              'Karachay-Cherkessia'),
    collection_site = factor(collection_site,
                             levels = c('Kabardino-Balkaria',
                                        'Karachay-Cherkessia'))
  )

# Order tips top-to-bottom
efaecalis_tip_order <- efaecalis_tip_coords |>
  arrange(desc(y)) |>
  mutate(tip_index = row_number()) |>
  select(label, tip_index)


# Plot the data -----------------------------------------------------------

efaecalis_tip_label_offset <- 0.025 *
  diff(range(efaecalis_edges$x_parent,
             efaecalis_edges$x_child, na.rm = TRUE))

efaecalis_tree_x_range <- diff(range(efaecalis_edges$x_parent,
                                     efaecalis_edges$x_child, na.rm = TRUE))

efaecalis_leader_x_end <- max(efaecalis_edges$x_child, na.rm = TRUE) +
  (0.12 * efaecalis_tree_x_range)

efaecalis_label_x <- efaecalis_leader_x_end +
  (0.02 * efaecalis_tree_x_range)

efaecalis_dot_start_offset <- 0.04 * efaecalis_tree_x_range

efaecalis_tree_left_x <- min(efaecalis_edges$x_parent,
                        efaecalis_edges$x_child, na.rm = TRUE) -
  (0.05 * efaecalis_tree_x_range)

# Tree segments
efaecalis_tree_plot <- ggplot() +
  geom_segment(data = efaecalis_edges_h,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.5) +
  geom_segment(data = efaecalis_edges_v,
               aes(x = x, y = y, xend = xend, yend = yend),
               linewidth = 0.5) +
  
  # Caucasus Region annotation
  geom_point(
    data = efaecalis_tip_coords |> left_join(efaecalis_tip_order, by = 'label'),
    aes(x = efaecalis_tree_left_x, y = y,
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
    values = efaecalis_site_cols
  ) +
  guides(
    shape = guide_legend(order = 1,
                         override.aes = list(size = 4)),
    color = guide_legend(order = 1)
  ) +
  
  ggnewscale::new_scale_color() +
  
  # Cluster legend
  geom_point(
    data = efaecalis_tip_coords |> left_join(efaecalis_tip_order, by = 'label'),
    aes(x = x + efaecalis_tip_label_offset * 0.3,
        y = y, color = cluster),
    size = 2
  ) +
  geom_segment(
    data = efaecalis_tip_coords |> left_join(efaecalis_tip_order, by = 'label'),
    aes(x = x + efaecalis_dot_start_offset,
        y = y,
        xend = efaecalis_leader_x_end,
        yend = y),
    linetype = 'dotted',
    linewidth = 0.4,
    color = 'grey70'
  ) +
  geom_text(
    data = efaecalis_tip_coords |> left_join(efaecalis_tip_order, by = 'label'),
    aes(x = efaecalis_label_x,
        y = y,
        label = label,
        color = cluster),
    hjust = 0,
    size = 3.5,
    show.legend = FALSE
  ) +
  scale_color_manual(
    name   = 'Clusters',
    values = efaecalis_cluster_cols
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
efaecalis_tree_plot


# Combine Tree and Phylogenetic Distance ----------------------------------

efaecalis_combined_tree_dist_plot <-
  efaecalis_tree_plot + efaecalis_hm +
  plot_layout(widths = c(1, 1.2), guides = 'collect') +
  plot_annotation(
    title = '(C)',
    subtitle = expression(bolditalic('Enterococcus faecalis')),
    theme = theme(
      plot.title = element_text(face = 'bold', size = 12),
      plot.subtitle = element_text(size = 12, hjust = 0.5),
      legend.title = element_text(face = 'bold', size = 12),
      legend.text  = element_text(size = 12)
    )
  )

# Check the plot
efaecalis_combined_tree_dist_plot

# Save the plot
ggsave(
  efaecalis_combined_tree_dist_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/efaecalis_combined_plot.png',
  height = 15,
  width = 32,
  units = 'cm',
  dpi = 300
)

# Save the plot
ggsave(
  efaecalis_combined_tree_dist_plot,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/efaecalis_combined_plot.tiff',
  height = 15,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine the Plots of Supplementary Bacteria -----------------------------

combined_strain_plot_supp <-
  (wrap_elements(full = lparacasei_combined_tree_dist_plot) /
     wrap_elements(full = lparakefiri_combined_tree_dist_plot) /
     wrap_elements(full = efaecalis_combined_tree_dist_plot))

# Check the plot
combined_strain_plot_supp

# Save the plot
ggsave(
  combined_strain_plot_supp,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/combined_strain_plot_supplementary.png',
  height = 40,
  width = 32,
  units = 'cm',
  dpi = 300
)

# Save the plot
ggsave(
  combined_strain_plot_supp,
  filename = 'plots/06_strain_profiling/6.1_strainphlan_tree_dist/combined_strain_plot_supplementary.tiff',
  height = 40,
  width = 32,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)
