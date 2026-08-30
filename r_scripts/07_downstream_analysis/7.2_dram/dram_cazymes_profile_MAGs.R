# DRAM MAGs Functional Profile
# ~/Documents/RStudio/mk_caucasus/scripts/07_downstream_analysis/7.2_dram
# 02/10/2025

# This is the Revised and Working Script ----------------------------------

# Load the necessary packages.
library(paletteer)
library(patchwork)
library(tidyverse)


# Load and Check Data -----------------------------------------------------

# Load data
dram_cazymes_data_full <- read_tsv('data/07_downstream_analysis/7.2_dram/hq_mags/distilled/combined_metabolism_summary.tsv', show_col_types = FALSE) |> 
  filter(gene_number > 0)

dram_cazymes_data <- read_tsv('data/07_downstream_analysis/7.2_dram/hq_mags/distilled/combined_metabolism_summary.tsv', show_col_types = FALSE) |> 
  select(-gene_description, -header, -subheader) |> 
  filter(gene_number > 0)

# Load gtdb data as well for species annotation
gtdb_data <- read_tsv('data/07_downstream_analysis/7.1_gtdbtk/hq_mags/gtdb_hqmags_with_qc_data.tsv', show_col_types = FALSE) |> 
  select(user_genome, species, sample_type, collection_site)

# Final data (combine dram and gtdb data)
final_cazymes_data <- dram_cazymes_data |>
  left_join(gtdb_data, by = c('bin_name' = 'user_genome'))

# Determine which module has the highest frequency (will be printed in the console)
final_cazymes_data |> 
  count(module, sort = TRUE)

# Subset glycoside hydrolases (gh)
cazymes_gh_data <- final_cazymes_data |> 
  filter(module == 'Glycoside Hydrolases')

# Determine which gene has the highest frequency (will be printed in the console)
cazymes_gh_data |> 
  count(gene_id, sort = TRUE)

# Subset glycosyl transferases (gt)
cazymes_gt_data <- final_cazymes_data |> 
  filter(module == 'Glycosyl Transferases')

# Determine which gene has the highest frequency (will be printed in the console)
cazymes_gt_data |> 
  count(gene_id, sort = TRUE)


# Glycoside Hydrolases ----------------------------------------------------

# Arrange bin_name and gene_id numerically and alphabetically
cazymes_gh_data <- cazymes_gh_data |> 
  # arrange bin_name
  mutate(
    sample_id = str_extract(bin_name, '^k\\d+_[a-z]'), # get something like "k1_g" from "k1_g_bin.1.permissive"
    k_num     = as.integer(str_extract(sample_id, '(?<=k)\\d+')), # extract the number after 'k'
    type      = str_extract(sample_id, '(?<=_)\\w+') # extract g / m (or whatever letter is after the underscore)
  ) |>
  # species frequency (unique bins per species)
  add_count(species, wt = !duplicated(bin_name), name = 'freq') |>
  # order: highest freq species first, ties alphabetically, then your bin order
  arrange(desc(freq), species, k_num, type, bin_name) |>
  mutate(
    bin_name = factor(bin_name, levels = rev(unique(bin_name)))
  ) |>
  # arrange gene_id
  mutate(
    gh_num = as.numeric(str_extract(gene_id, "\\d+"))
  ) |>
  arrange(gh_num) |>
  mutate(
    gene_id = factor(gene_id, levels = rev(unique(gene_id)))
  )

# Adding color to species annotation (match the plotting order)
gh_species_levels <- cazymes_gh_data |>
  distinct(species, bin_name) |>
  count(species, name = 'freq') |>
  arrange(desc(freq), species) |>
  pull(species)

gh_species_color  <- paletteer_d('ggthemes::Tableau_10') |> 
  as.character()
gh_species_color  <- gh_species_color[seq_along(gh_species_levels)]

cazymes_gh_data <- cazymes_gh_data |>
  mutate(species = factor(species, levels = gh_species_levels))

# Plot the data
gh_bubble_plot <- ggplot(
  cazymes_gh_data,
  aes(x = bin_name, y = gene_id, size = gene_number, color = species)
) +
  geom_point(alpha = 0.8) +
  # Manual size scale: 1, 3, 5, 7, 9
  scale_size_continuous(
    name = 'Gene Count',
    breaks = c(1, 3, 5, 7, 9),
    range  = c(1, 10), # adjust bubble sizes if needed
    limits = c(1, 9),
    guide  = guide_legend(order = 2) # Gene Count legend second
  ) +
  scale_color_manual(
    name   = 'Species',
    values = gh_species_color,
    guide  = guide_legend(
      order       = 1, # Species legend on top
      title.theme = element_text(face = 'bold', size = 10),
      label.theme = element_text(face = 'italic', size = 10),
      override.aes = list(size = 5) # make legend circles bigger
    )
  ) +
  scale_y_discrete(
    limits = sort(unique(cazymes_gh_data$gene_id), decreasing = TRUE),
    expand = expansion(mult = c(0.022, 0.022)) # add vertical space
  ) +
  scale_x_discrete(position = 'top',
    expand = expansion(mult = c(0.022, 0.022))
  ) +
  coord_flip() +
  labs(
    title = '(A)',
    subtitle = 'Glycoside Hydrolase (n = 759, gene_IDs = 56)',
    x = 'High Quality MAGs',
    y = 'Glycoside Hydrolase Gene ID',
    size = 'Gene Number'
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(face = 'bold', size = 12),
    plot.subtitle    = element_text(face = 'bold', size = 12, hjust = 0.5),
    axis.title       = element_text(face = 'bold', size = 12),
    axis.title.x     = element_text(margin = margin(t = 10)),
    axis.text.x      = element_text(size = 10, angle = 45, hjust = 1),
    # axis.title.y     = element_text(margin = margin(r = 10)),
    axis.title.y     = element_blank(),
    axis.text.y      = element_blank(),
    #axis.ticks.y     = element_lin(),
    # legend.title   = element_text(face = 'bold', size = 10),
    legend.position  = 'none', # removes the legend so there's only one legend in GT
    panel.grid.major = element_line(linetype = 'dashed', color = 'grey90', linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Check the plot
gh_bubble_plot

# Save the plot
ggsave(
  gh_bubble_plot,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazymes_gh.tiff',
  height = 30,
  width = 40,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Glycosyl Transferases ---------------------------------------------------

# Arrange bin_name and gene_id numerically and alphabetically
cazymes_gt_data <- cazymes_gt_data |> 
  # arrange bin_name
  mutate(
    sample_id = str_extract(bin_name, '^k\\d+_[a-z]'), # get something like "k1_g" from "k1_g_bin.1.permissive"
    k_num     = as.integer(str_extract(sample_id, '(?<=k)\\d+')), # extract the number after 'k'
    type      = str_extract(sample_id, '(?<=_)\\w+') # extract g / m (or whatever letter is after the underscore)
  ) |>
  # species frequency (unique bins per species)
  add_count(species, wt = !duplicated(bin_name), name = 'freq') |>
  # order: highest freq species first, ties alphabetically, then your bin order
  arrange(desc(freq), species, k_num, type, bin_name) |>
  mutate(
    bin_name = factor(bin_name, levels = rev(unique(bin_name)))
  ) |>
  # arrange gene_id
  mutate(
    gt_num = as.numeric(str_extract(gene_id, "\\d+"))
  ) |>
  arrange(gt_num) |>
  mutate(
    gene_id = factor(gene_id, levels = rev(unique(gene_id)))
  )

# Adding color to species annotation (match the plotting order)
gt_species_levels <- cazymes_gt_data |>
  distinct(species, bin_name) |>
  count(species, name = 'freq') |>
  arrange(desc(freq), species) |>
  pull(species)

gt_species_color  <- paletteer_d('ggthemes::Tableau_10') |> 
  as.character()
gt_species_color  <- gt_species_color[seq_along(gt_species_levels)]

cazymes_gt_data <- cazymes_gt_data |>
  mutate(species = factor(species, levels = gt_species_levels))

# Plot the data
gt_bubble_plot <- ggplot(
  cazymes_gt_data,
  aes(x = bin_name, y = gene_id, size = gene_number, color = species)
) +
  geom_point(alpha = 0.8) +
  # Manual size scale: 1, 3, 5, 7, 9
  scale_size_continuous(
    name = 'Gene Count',
    breaks = c(2, 4, 6, 8, 10),
    range  = c(1, 10), # adjust bubble sizes if needed
    limits = c(1, 10),
    guide  = guide_legend(order = 2) # Gene Count legend second
  ) +
  scale_color_manual(
    name   = 'Species',
    values = gt_species_color,
    guide  = guide_legend(
      order       = 1, # Species legend on top
      title.theme = element_text(face = 'bold', size = 12),
      label.theme = element_text(face = 'italic', size = 12),
      override.aes = list(size = 5) # make legend circles bigger
    )
  ) +
  scale_y_discrete(
    limits = sort(unique(cazymes_gt_data$gene_id), decreasing = TRUE),
    expand = expansion(mult = c(0.04, 0.04)) # add vertical space
  ) +
  scale_x_discrete(expand = expansion(mult = c(0.022, 0.022))) +
  coord_flip() +
  labs(
    title = '(B)',
    subtitle = 'Glycosyl Transferases (n = 289, gene_IDs = 20)',
    x = 'High Quality MAGs',
    y = 'Glycosyl Transferases Gene ID',
    size = 'Gene Number'
  ) +
  theme_bw() +
  theme(
    plot.title       = element_text(face = 'bold', size = 12),
    plot.subtitle    = element_text(face = 'bold', size = 12, hjust = 0.5),
    axis.title       = element_text(face = 'bold', size = 12),
    axis.title.x     = element_text(margin = margin(t = 10)),
    # axis.title.y   = element_text(margin = margin(r = 10)),
    axis.title.y     = element_blank(), # removes the y-axis title so there's only one legend in GT
    axis.text.x      = element_text(size = 10, angle = 45, hjust = 1),
    axis.text.y      = element_blank(), # removes the y-axis text so there's only one legend in GT
    legend.title     = element_text(face = 'bold', size = 12),
    legend.text      = element_text(size = 12),
    panel.grid.major = element_line(linetype = 'dashed', color = 'grey90', linewidth = 0.3),
    panel.grid.minor = element_blank()
  )

# Check the plot
gt_bubble_plot

# Save the plot
ggsave(
  gt_bubble_plot,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazymes_gt.tiff',
  height = 30,
  width = 30,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)


# Combine GH and GT Plot --------------------------------------------------

library(gridExtra)

cazymes_combined_plot <- grid.arrange(
  gh_bubble_plot, 
  gt_bubble_plot, 
  ncol = 2,
  widths = c(1.6, 1)
  )

# Save the plot
ggsave(
  cazymes_combined_plot,
  filename = 'plots/07_downstream_analysis/7.2_dram/cazymes_combined_plot.tiff',
  height = 35,
  width = 50,
  units = 'cm',
  dpi = 600,
  compression = 'lzw'
)