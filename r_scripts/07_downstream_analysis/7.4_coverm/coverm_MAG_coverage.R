# CoverM MAG Coverage
# ~/Documents/RStudio/mk_caucasus/scripts/07_downstream_analysis/7.4_coverm
# 16/01/2025

# Load the necessary packages.
library(tidyverse)


# Load and Check Data -----------------------------------------------------

coverm_data <- read_tsv('data/07_downstream_analysis/7.4_coverm/coverm_mag_coverage_summary.tsv', show_col_types = FALSE)


# Setting of Threshold for Prevalence -------------------------------------

# Choose a presence threshold for covered_fraction
# Common choices:
#   0      = any coverage at all
#   0.1    = more conservative "present"
#   0.5    = strong presence
cf_thresh <- 0.1


# Making a Summary File ---------------------------------------------------

coverm_summary <- coverm_data |>
  rowwise() |>
  mutate(
    # max coverage based on mean depth
    coverage_max = max(c_across(ends_with('_mean')), na.rm = TRUE),
    
    # prevalence based on covered fraction threshold
    prevalence = sum(c_across(ends_with('_covered_fraction')) >= cf_thresh, na.rm = TRUE),
    
    # optional: max covered fraction across samples (nice for an extra ring/track)
    covered_fraction_max = max(c_across(ends_with('_covered_fraction')), na.rm = TRUE)
  ) |>
  ungroup() |>
  select(
    bin = genome,
    coverage_max,
    prevalence,
    covered_fraction_max
  )

# Save data (to be merged in phylophlan metadata)
write_tsv(coverm_summary, 'data/07_downstream_analysis/7.4_coverm/coverm_summary.tsv')
