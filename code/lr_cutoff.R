# Load Required Libraries
suppressMessages(suppressWarnings({
  library(dplyr)
  library(furrr)
  library(data.table)
}))

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
slurm_job_id <-  as.character(args[1])

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

output_dir <- file.path("output", paste0("simulation_", slurm_job_id))

# load individual lrs and assemble
t<-getwd()
setwd(paste0(t,"/temp"))
tabs<-list.files(pattern = "sim_processed_genotypes", recursive = TRUE)
final_results<-do.call(rbind, lapply(tabs, read.csv))

# Define Populations
population_labels <- tibble(
  population = factor(
    c("AfAm", "Cauc", "Hispanic", "Asian"),
    levels = c("AfAm", "Cauc", "Hispanic", "Asian")
  ),
  label = c("African American", "Caucasian", "Hispanic", "Asian"),
  population_number = c(1:4)
)

#Assemble individual genotypes
final_results<-left_join(final_results,as.data.frame(population_labels))
setwd(output_dir)
fwrite(final_results,"sim_processed_genotypes.csv")
rm(final_results)

# load combined lrs and assemble
# Assemble summary genotypes
setwd(paste0(t,"/temp"))
tabs<-list.files(pattern = "sim_combined_genotypes", recursive = TRUE)
combined_lrs<-do.call(rbind, lapply(tabs, read.csv))

#To load if pre-assembled:
#combined_lrs = fread(paste0("output/simulation_",slurm_job_id,"/sim_summary_genotypes.csv"))

# Function to calculate cut-off values for 1%, 0.1%, and 0.01% FPR
calculate_cutoffs <- function(input_df, fp_rates) {
  cutoffs <- input_df |>
    filter(relationship_known == "unrelated") |>
    group_by(population, loci_set) |>
    summarize(
      fixed_cutoff = 1.00,
      cutoff_1 = quantile(LR, probs = 1 - fp_rates[1] / 100, na.rm = TRUE),
      cutoff_0_1 = quantile(LR, probs = 1 - fp_rates[2] / 100, na.rm = TRUE),
      cutoff_0_01 = quantile(LR, probs = 1 - fp_rates[3] / 100, na.rm = TRUE),
      n_unrelated = n(),
      .groups = 'drop'
    )
  return(cutoffs)
}

calculate_proportions_exceeding_cutoffs <- function(input_df, cutoffs) {
  df_with_cutoffs <- left_join(input_df, cutoffs, by = c("population", "loci_set"))
  df_with_cutoffs <- df_with_cutoffs |>
    mutate(
      exceeds_fixed_cutoff = LR > fixed_cutoff,
      exceeds_cutoff_1 = LR > cutoff_1,
      exceeds_cutoff_0_1 = LR > cutoff_0_1,
      exceeds_cutoff_0_01 = LR > cutoff_0_01
    )
  proportions_exceeding <- df_with_cutoffs |>
    group_by(population, relationship_tested, loci_set) |>
    summarize(
      proportion_exceeding_fixed = sum(exceeds_fixed_cutoff, na.rm = TRUE) / n(),
      proportion_exceeding_1 = sum(exceeds_cutoff_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_1 = sum(exceeds_cutoff_0_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_01 = sum(exceeds_cutoff_0_01, na.rm = TRUE) / n(),
      n_related = n(),
      .groups = 'drop'
    ) |>
    filter(relationship_tested != "unrelated")
  return(proportions_exceeding)
}

# Calculate and save cutoffs
cutoffs <- calculate_cutoffs(combined_lrs, c(1, 0.1, 0.01))
fwrite(cutoffs, file.path(output_dir, "sim_cutoffs.csv"))
    
proportions_exceeding_cutoffs <- calculate_proportions_exceeding_cutoffs(combined_lrs, cutoffs)
fwrite(proportions_exceeding_cutoffs, file.path(output_dir, "sim_proportions_exceeding_cutoffs.csv"))
log_message("Cutoffs applied.")
