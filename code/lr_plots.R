# Load Required Libraries
suppressMessages(suppressWarnings({
  library(ggplot2)
  library(dplyr)
  library(tibble)
  library(furrr)
  library(tidyr)
  library(data.table)
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

file_type<-"png" #can be pdf or png, pdf is faster
t<-getwd()
output_dir <- file.path("output", paste0("simulation_", slurm_job_id))

# load combined_lrs
combined_lrs = fread(paste0("output/simulation_",slurm_job_id,"/sim_summary_genotypes.csv"))
combined_lrs$LR<-as.numeric(combined_lrs$LR)
combined_lrs$relationship_tested <- factor(combined_lrs$relationship_tested, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
combined_lrs$population<-as.factor(combined_lrs$population)
# load proportions_exceeding_cutoffs
proportions_exceeding_cutoffs = fread(paste0("output/simulation_",slurm_job_id,"/sim_proportions_exceeding_cutoffs.csv"))

summary_stats <- combined_lrs[,    list(
  mean_LR = mean(LR),
  lower_95 = quantile(LR, 0.025),
  upper_95 = quantile(LR, 0.975)#,
),
by=c("relationship_tested","population","loci_set")]

plot_and_save_results <- function(combined_lrs,file_type) {
  #New option added to print images to pdf files if you don't need high resolution pngs. Writing to pdf will be faster in development phases.

  log_message("Starting to plot LR distributions...")
  print(head(combined_lrs))
  # Ensure factor levels are set correctly for plotting

  plot1<-ggplot(combined_lrs, aes(x = relationship_tested, y = LR, fill = population, color = population)) +
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

    if(!file_type=="png"){
      pdf(paste(output_dir, "/sim_log_lr_panel_plot.pdf", sep = ""))
      invisible(print(plot1))
      dev.off()
    } else {
      ggsave(plot=plot1,path=paste0(output_dir, "/"),filename="sim_log_lr_panel_plot.png", width = 12, height = 8)
    }

    log_message("LR distributions plot saved.")


    plot2<-ggplot(summary_stats, aes(x = loci_set, y = mean_LR, group = population, color = population)) +
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

    if(!file_type=="png"){
      pdf(paste(output_dir, "/sim_line_chart_lr.pdf", sep = ""))
      invisible(print(plot2))
      dev.off()
    } else {
      ggsave(plot=plot2,path=paste0(output_dir, "/"),filename="sim_line_chart_lr.png", width = 14, height = 10)
    }

    log_message("Mean LR plot saved.")
}

plot_relationship_comparisons <- function(combined_lrs,file_type) {
  log_message("Starting to plot relationship comparisons...")

  plot3<-ggplot(combined_lrs, aes(x = relationship_tested, y = LR, fill = population)) +
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


  if(!file_type=="png"){
    pdf(paste(output_dir, "/relationship_comparisons_plot.pdf", sep = ""))
    invisible(print(plot3))
    dev.off()
  } else {
  #  setwd(output_dir)
    ggsave(plot=plot3,path=paste0(output_dir, "/"),filename="relationship_comparisons_plot.png", width = 16, height = 12)
  }

  log_message("Relationship comparisons plot saved.")
}

plot_proportions_exceeding_cutoffs <- function(proportions_exceeding_cutoffs,file_type) {
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

  plot4<-ggplot(proportions_long, aes(x = relationship_tested, y = Proportion, fill = population, color = population)) +
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

  if(!file_type=="png"){
    pdf(paste(output_dir, "/relationship_comparisons_plot.pdf", sep = ""))
    print(plot4)
    dev.off()
  } else {
    ggsave(plot=plot4,path=paste0(output_dir, "/"),filename="sim_proportions_exceeding_cutoffs_combined.png", width = 12, height = 8)
  }

  log_message("Proportions exceeding cutoffs plot saved.")
}

plot_and_save_results(combined_lrs,file_type)

plot_relationship_comparisons(combined_lrs,file_type)

plot_proportions_exceeding_cutoffs(proportions_exceeding_cutoffs,file_type)

closeAllConnections()

