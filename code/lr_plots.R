
# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(furrr)
  library(data.table)
  library(ggplot2)
  library(future)
  library(parallel)
}))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
slurm_job_id <- as.character(args[1])

output_dir <- file.path("output", paste0("simulation_", slurm_job_id))

# load combined_lrs
combined_lrs = fread(paste0("output/simulation_",slurm_job_id,"/sim_summary_genotypes.csv"))

# load proportions_exceeding_cutoffs
proportions_exceeding_cutoffs = fread(paste0("output/simulation_",slurm_job_id,"/sim_proportions_exceeding_cutoffs.csv"))

plot_and_save_results <- function(combined_lrs) {
  log_message("Starting to plot LR distributions...")
  
  # Ensure factor levels are set correctly for plotting
  combined_lrs$relationship_tested <- factor(combined_lrs$relationship_tested, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
  
  ggplot(combined_lrs, aes(x = relationship_tested, y = LR, fill = population, color = population)) +
    geom_boxplot(position = position_dodge(width = 0.9)) +
    facet_wrap(~ loci_set, scales = "fixed") +
    labs(
      title = "LR Distributions Across Populations and Relationship Types",
      x = "Relationship Tested",
      y = "LR",
      fill = "Population",
      color = "Population"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_log10()
  
  ggsave(file.path(output_dir, "sim_log_lr_panel_plot.png"), width = 12, height = 8)
  
  log_message("LR distributions plot saved.")
  
  summary_stats <- combined_lrs |>
    group_by(relationship_tested, population, loci_set) |>
    summarize(
      mean_LR = mean(LR),
      lower_95 = quantile(LR, 0.025),
      upper_95 = quantile(LR, 0.975),
      .groups = 'drop'
    )
  
  ggplot(summary_stats, aes(x = loci_set, y = mean_LR, group = population, color = population)) +
    geom_line(size = 1) +
    facet_wrap(~ relationship_tested, scales = "free_y", ncol = 2) +
    scale_y_log10() +
    labs(
      title = "Mean LR Across Populations and Relationship Types",
      x = "Loci Set",
      y = "Combined LR",
      color = "Population"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1),
      legend.position = "bottom"
    )
  
  ggsave(file.path(output_dir, "sim_line_chart_lr.png"), width = 14, height = 10)
  
  log_message("Mean LR plot saved.")
}

plot_relationship_comparisons <- function(combined_lrs) {
  log_message("Starting to plot relationship comparisons...")
  
  ggplot(combined_lrs, aes(x = relationship_tested, y = LR, fill = population)) +
    geom_boxplot(position = position_dodge(width = 0.9)) +
    facet_grid(relationship_known ~ loci_set) +
    labs(
      title = "LR Distributions for Relationship Comparisons",
      x = "Tested Relationship Type",
      y = "LR",
      fill = "Population"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_log10()
  
  ggsave(file.path(output_dir, "relationship_comparisons_plot.png"), width = 16, height = 12)
  
  log_message("Relationship comparisons plot saved.")
}

plot_proportions_exceeding_cutoffs <- function(proportions_exceeding_cutoffs) {
  log_message("Starting to plot proportions exceeding cutoffs...")
  
  # Ensure factor levels are set correctly for plotting
  proportions_exceeding_cutoffs$relationship_tested <- factor(proportions_exceeding_cutoffs$relationship_tested, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
  proportions_exceeding_cutoffs$population <- factor(proportions_exceeding_cutoffs$population, levels = c("AfAm", "Cauc", "Hispanic", "Asian"))
  
  proportions_long <- proportions_exceeding_cutoffs |>
    pivot_longer(cols = starts_with("proportion_exceeding"),
                 names_to = "Cutoff_Type", values_to = "Proportion",
                 names_prefix = "proportion_exceeding_")
  
  proportions_long$Cutoff_Type <- factor(proportions_long$Cutoff_Type, levels = c("fixed", "1", "0_1", "0_01"),
                                         labels = c("Fixed Cutoff (1.00)", "1% FPR", "0.1% FPR", "0.01% FPR"))
  
  ggplot(proportions_long, aes(x = relationship_tested, y = Proportion, fill = population, color = population)) +
    geom_bar(stat = "identity", position = position_dodge(width = 0.9)) +
    facet_wrap(~ Cutoff_Type + loci_set, scales = "fixed") +
    labs(
      title = "Proportions Exceeding Likelihood Cut-offs Across Relationship Types and Loci Sets",
      x = "Relationship Tested",
      y = "Proportion Exceeding Cut-off",
      fill = "Population",
      color = "Population"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    )
  
  ggsave(file.path(output_dir, "sim_proportions_exceeding_cutoffs_combined.png"), width = 12, height = 8)
  
  log_message("Proportions exceeding cutoffs plot saved.")
}


plot_and_save_results(combined_lrs)
plot_relationship_comparisons(combined_lrs)
plot_proportions_exceeding_cutoffs(proportions_exceeding_cutoffs)


