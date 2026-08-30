# PhyloPhlAn Phylogenetic Tree
# ~/Documents/RStudio/mk_caucasus/scripts/07_downstream_analysis/7.5_phylophlan
# 14/01/2026

# Load the necessary packages.

library(ape)
library(ggnewscale)
library(ggtree)
library(grid)
library(paletteer)
library(scales)
library(tidyverse)


# Prepare and Clean the data ----------------------------------------------

# Load phylophlan metadata
phylophlan_metadata <- read_tsv('data/07_downstream_analysis/7.5_phylophlan/hq_mags/phylophlan_hq_mags_metadata.tsv', show_col_types = FALSE)

# Load coverm metadata
coverm_metadata <- read_tsv('data/07_downstream_analysis/7.6_coverm/coverm_summary.tsv', show_col_types = FALSE)

# Combine: keep ALL rows from phylophlan, add matching coverm columns
combined_metadata <- phylophlan_metadata |>
  left_join(coverm_metadata, by = c('user_genome' = 'bin'))

# Sanity checks (highly recommended)
# 1) Any phylophlan genomes that didn't find a match in coverm?
unmatched <- combined_metadata |>
  filter(if_any(-user_genome, ~ FALSE)) # placeholder to keep pipe valid

setdiff(phylophlan_metadata$user_genome, coverm_metadata$bin)

# 2) Any duplicates in keys? (should be none)
phylophlan_metadata |> count(user_genome) |> filter(n > 1)
coverm_metadata |> count(bin) |> filter(n > 1)

# Save the final metadata
write_tsv(combined_metadata, 'data/07_downstream_analysis/7.5_phylophlan/hq_mags/phylophlan_hq_mags_final_metadata.tsv')


# Load the Final Metadata and the Tree ------------------------------------

# Delete the existing objects in rstudio

metadata <- read_tsv('data/07_downstream_analysis/7.5_phylophlan/hq_mags/phylophlan_hq_mags_final_metadata.tsv', show_col_types = FALSE)

# Load phylophlan tree.
phylotree <- read.tree('data/07_downstream_analysis/7.5_phylophlan/hq_mags/output_tol/input_genomes.tre.treefile')


# Add Necessary Annotations -----------------------------------------------

# 1) Base tree
p_base <- ggtree(phylotree, layout = 'circular')

# 2) Tip coordinates from the plotted tree
tips <- p_base$data |>
  filter(isTip) |>
  select(label, x, y)

# 3) Join your metadata
ann <- metadata |>
  transmute(
    label = user_genome,
    species,
    completeness,
    contamination,
    coverage_max,
    prevalence,
    closest_genome_ani
  ) |>
  left_join(tips, by = 'label') |>
  drop_na(x, y) |>
  mutate(
    closest_genome_ani = as.numeric(closest_genome_ani),
    species_classification = if_else(
      is.na(closest_genome_ani) | closest_genome_ani < 95,
      'Putative Novel Species (<95% ANI)',
      'Known Species (>95% ANI)'
    ),
    # force both legend levels to exist
    species_classification = factor(
      species_classification,
      levels = c('Known Species (>95% ANI)', 'Putative Novel Species (<95% ANI)')
    )
  )

# 4) Define ring positions (x positions)
max_x <- max(ann$x)

ann <- ann |>
  mutate(
    ring_comp_x  = max_x + 0.30,  # completeness (inner)
    ring_cont_x  = max_x + 0.45,  # contamination
    ring_cov_x   = max_x + 0.60,  # coverage_max
    ring_prev_x  = max_x + 0.75,  # prevalence
    ring_ani_x   = max_x + 0.90   # species classification (outermost)
  )

# Dummy legend data (forces BOTH brown + gold to show, even if absent)
legend_df <- tibble(
  x = ann$ring_ani_x[1],
  y = ann$y[1],
  species_classification = factor(
    c('Known Species (>95% ANI)', 'Putative Novel Species (<95% ANI)'),
    levels = levels(ann$species_classification)
  )
)


# Plotting the Tree with Annotations --------------------------------------

# 5) Expand plot limits to make room for all rings
annotated_phylotree_plot <- p_base +
  xlim_tree(max_x + 1.05) +
  
  # Species dots at tips
  geom_point(
    data = ann,
    aes(x = x, y = y, color = species),
    size = 3
  ) +
  scale_color_paletteer_d(
    palette = 'ggthemes::Classic_10',
    guide = guide_legend(
      title = 'Species',
      ncol = 2,
      order = 1,
      label.theme = element_text(face = 'italic')
    )
  ) +
  
  # Ring 1: Completeness (fill gradient; continuous bar)
  geom_tile(
    data = ann,
    aes(x = ring_comp_x, y = y, fill = completeness),
    width = 0.12,
    height = 1
  ) +
  scale_fill_gradient(
    name = 'Completeness (%)',
    limits = c(90, 100),
    low = 'white',
    high = '#008000',
    oob = scales::squish,
    guide = guide_colorbar(
      order = 2,
      direction = 'horizontal',
      title.position = 'top',
      label.position = 'bottom',
      barwidth  = grid::unit(6, 'cm'),
      barheight = grid::unit(0.35, 'cm')
    )
  ) +
  
  # Ring 2: Contamination (fill; continuous bar)
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = ann,
    aes(x = ring_cont_x, y = y, fill = contamination),
    width = 0.12,
    height = 1
  ) +
  scale_fill_gradient(
    name = 'Contamination (%)',
    limits = c(0, 5),
    low = 'white',
    high = 'red',
    oob = scales::squish,
    guide = guide_colorbar(
      order = 3,
      direction = 'horizontal',
      title.position = 'top',
      label.position = 'bottom',
      barwidth  = grid::unit(5.5, 'cm'),
      barheight = grid::unit(0.30, 'cm')
    )
  ) +
  
  # Ring 3: Coverage_max (fill; continuous bar)
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = ann,
    aes(x = ring_cov_x, y = y, fill = coverage_max),
    width = 0.12,
    height = 1
  ) +
  scale_fill_gradient(
    name = 'Coverage max',
    limits = c(0, 120),
    low = 'white',
    high = '#008080',
    oob = scales::squish,
    guide = guide_colorbar(
      order = 4,
      direction = 'horizontal',
      title.position = 'top',
      label.position = 'bottom',
      barwidth  = grid::unit(5.5, 'cm'),
      barheight = grid::unit(0.30, 'cm')
    )
  ) +
  
  # Ring 4: Prevalence (fill; continuous bar)
  ggnewscale::new_scale_fill() +
  geom_tile(
    data = ann,
    aes(x = ring_prev_x, y = y, fill = prevalence),
    width = 0.12,
    height = 1
  ) +
  scale_fill_gradient(
    name = 'Prevalence',
    limits = c(0, 21),
    low = 'white',
    high = '#800080',
    oob = scales::squish,
    guide = guide_colorbar(
      order = 5,
      direction = 'horizontal',
      title.position = 'top',
      label.position = 'bottom',
      barwidth  = grid::unit(5.5, 'cm'),
      barheight = grid::unit(0.30, 'cm')
    )
  ) +
  
  # Ring 5: Species Classification (BROWN/GOLD)
  ggnewscale::new_scale_fill() +
  
  # Actual outer ring as tiles (no outline)
  geom_tile(
    data = ann,
    aes(x = ring_ani_x, y = y, fill = species_classification),
    width = 0.12,
    height = 1
  ) +
  
  # Dummy layer INSIDE plot range so it won't shrink the tree
  # alpha=0 hides it on the plot, but legend keys will be forced visible via override.aes
  geom_tile(
    data = legend_df,
    aes(x = x, y = y, fill = species_classification),
    width = 0.12,
    height = 1,
    alpha = 0,
    inherit.aes = FALSE,
    show.legend = TRUE
  ) +
  
  scale_fill_manual(
    name = 'Species Classification',
    values = c(
      'Known Species (>95% ANI)' = '#E5BA41',
      'Putative Novel Species (<95% ANI)' = '#7B542F'
    ),
    breaks = c('Known Species (>95% ANI)', 'Putative Novel Species (<95% ANI)'),
    drop = FALSE,
    guide = guide_legend(
      order = 6,
      direction = 'vertical',
      ncol = 1,
      byrow = TRUE,
      title.position = 'top',
      # THIS is the key: force legend tiles to be opaque squares
      override.aes = list(alpha = 1)
    )
  ) +
  
  labs(
    title = '(C)'
  ) +
  theme_void() +
  theme(
    legend.position = 'right',
    legend.box = 'vertical',
    legend.title = element_text(face = 'bold', size = 12),
    legend.text = element_text(size = 12),
    legend.spacing.y = grid::unit(0.35, 'cm'),
    legend.key.width  = grid::unit(0.7, 'cm'),
    legend.key.height = grid::unit(0.35, 'cm'),
    plot.title = element_text(face = 'bold', size = 12, margin = margin(b = -60, l = 60)),
    plot.margin = margin(t = 10, r = 50, b = 0, l = 0)
  )

# Check the plot
annotated_phylotree_plot

# Save the plot
ggsave(
  annotated_phylotree_plot,
  filename = 'plots/07_downstream_analysis/7.5_phylophlan/annotated_phylotree_hq_mags.tiff',
  height = 20,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine CheckM2, GTDB-Tk, and PhyloPhlAn Outputs ------------------------

library(gridExtra)

# Option 1
layout_a <- rbind(
  c(1, 2), # top plot spans 2 columns
  c(3, 3) # bottom plot spans 2 columns
)

mags_combined <- grid.arrange(
  gtdb_lollipop_plot,
  gtdb_species_count_a,
  annotated_phylotree_plot,
  layout_matrix = layout_a,
  heights = c(1, 1.5),
  widths = c(1.2, 1)
)

# Save the plot
ggsave(
  mags_combined,
  filename = 'plots/07_downstream_analysis/7.5_phylophlan/mags_combined_plot.tiff',
  height = 35,
  width = 40,
  units = 'cm',
  bg = 'white',
  dpi = 600,
  compression = 'lzw'
)