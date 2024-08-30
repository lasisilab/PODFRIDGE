# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(furrr)
  library(data.table)
  library(future)
  library(parallel)
}))

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
slurm_job_id <-  as.character(args[1])

output_dir <- file.path("output", paste0("simulation_", slurm_job_id))

# load combined lrs
combined_lrs = fread(paste0("output/simulation_",slurm_job_id,"/sim_summary_genotypes.csv"))
summary(combined_lrs)
combined_lrs$LR<-as.numeric(combined_lrs$LR) ##

# Function to calculate cut-off values for 1%, 0.1%, and 0.01% FPR
calculate_cutoffs <- function(input_df, fp_rates) {
  cutoffs <- input_df[relationship_known == "unrelated",
                                  list(
      fixed_cutoff = 1.00,
      cutoff_1 = quantile(LR, probs = 1 - fp_rates[1] / 100, na.rm = TRUE),
      cutoff_0_1 = quantile(LR, probs = 1 - fp_rates[2] / 100, na.rm = TRUE),
      cutoff_0_01 = quantile(LR, probs = 1 - fp_rates[3] / 100, na.rm = TRUE),
      n_unrelated = .N),
    by=c("population", "loci_set")]
  return(cutoffs)
}

calculate_proportions_exceeding_cutoffs <- function(input_df, cutoffs) {
  df_with_cutoffs <- left_join(input_df, cutoffs, by = c("population", "loci_set"))
  df_with_cutoffs <- df_with_cutoffs [,list(
      population=population,
      relationship_tested = relationship_tested,
      loci_set = loci_set,
      exceeds_fixed_cutoff = LR > fixed_cutoff,
      exceeds_cutoff_1 = LR > cutoff_1,
      exceeds_cutoff_0_1 = LR > cutoff_0_1,
      exceeds_cutoff_0_01 = LR > cutoff_0_01
    )]

  df_no_cutoffs<-df_with_cutoffs[!df_with_cutoffs$relationship_tested == "unrelated",]
  proportions_exceeding <- df_no_cutoffs[,    list(
      proportion_exceeding_fixed = sum(exceeds_fixed_cutoff) / .N,
      proportion_exceeding_1 = sum(exceeds_cutoff_1, na.rm = TRUE) / .N,
      proportion_exceeding_0_1 = sum(exceeds_cutoff_0_1, na.rm = TRUE) / .N,
      proportion_exceeding_0_01 = sum(exceeds_cutoff_0_01, na.rm = TRUE) / .N,
      n_related = .N),
      by=c("population","relationship_tested","loci_set")]
  return(proportions_exceeding)
}

# Calculate and save cutoffs
cutoffs <- calculate_cutoffs(combined_lrs, c(1, 0.1, 0.01))
fwrite(cutoffs, file.path(output_dir, "sim_cutoffs.csv"))

proportions_exceeding_cutoffs <- calculate_proportions_exceeding_cutoffs(combined_lrs, cutoffs)
fwrite(proportions_exceeding_cutoffs, file.path(output_dir, "sim_proportions_exceeding_cutoffs.csv"))
