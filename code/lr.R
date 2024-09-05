# Load Required Libraries
library(tidyverse)
library(furrr)
library(data.table)
library(future)
library(parallel)

# Set up cluster

# For on one machine if not using cluster
if(!exists("use_remote_cluster")){ #Add this parameter to run script
  use_remote_cluster<-0
} 
if(use_remote_cluster==0){
  if( .Platform$OS.type == "windows" ){
    cl <- makeCluster(availableCores())  #specify how many cores to use
  } else { # use the fork cluster type on linux because its faster- not available for windows
    cl <- makeCluster(availableCores(),type="FORK")  #specify how many cores to use
  } 
  plan(cluster, workers=cl)
  # Ensure the cluster is stopped when the script exits
  on.exit(parallel::stopCluster(cl))
} else { #if sending to remote cluster(s)
  plan(cluster, workers = clustername1) #use URL if using an online cluster, multiple clusters can also be specified here  
} 

# Ensure the cluster is stopped when the script exits
on.exit(parallel::stopCluster(cl))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
slurm_job_id <-  as.character(args[1])

# Create output folder with SLURM job ID
output_dir <- file.path("output", paste0("simulation_", slurm_job_id))
dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "sim_processed_genotypes.csv")
summary_output_file <- file.path(output_dir, "sim_summary_genotypes.csv")
timing_log_file <- file.path(output_dir, "timing_log.csv")

# Log the start of the process
log_message("Starting simulation setup and processing...")

# Load Allele Frequencies Data
log_message("Loading allele frequencies data...")
allele_freq_time <- system.time({
  df_allelefreq <- fread("data/df_allelefreq_combined.csv")
  df_allelefreq <- df_allelefreq[population != "all"] # Filter out "all" population
  df_allelefreq[, allele := as.character(allele)]
})
log_message(paste("Loaded allele frequencies data in", allele_freq_time["elapsed"], "seconds."))

# Extract unique loci
log_message("Extracting unique loci...")
loci_list <- unique(df_allelefreq$marker)


# Load Core Loci Data
log_message("Loading core loci data...")
core_loci_time <- system.time({
  core_loci <- fread("data/core_CODIS_loci.csv")
  columns <- c("core_13", "identifiler_15", "expanded_20", "supplementary")
  loci_lists <- lapply(columns, function(col) {
    core_loci |> 
      filter(get(col) == 1) |> 
      pull(locus)
  })
  names(loci_lists) <- columns
  loci_lists$autosomal_29 <- loci_list
})
log_message(paste("Loaded core loci data in", core_loci_time["elapsed"], "seconds."))

# Load Individuals Genotypes Data
input_dir <- file.path("data", "sims", paste0("simulation_", job_id))
individuals_genotypes <- fread(paste0(input_dir,"/processed_genotypes.csv"))

# Define Kinship Matrix
kinship_matrix <- data.table(
  relationship_type = factor(c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins")),
  k0 = c(0, 1/4, 1/2, 7/8, 15/16),
  k1 = c(1, 1/2, 1/2, 1/8, 1/16),
  k2 = c(0, 1/4, 0, 0, 0)
)

# Define Populations
population_labels <- data.table(
  population = factor(c("AfAm", "Cauc", "Hispanic", "Asian")),
  label = c("African American", "Caucasian", "Hispanic", "Asian")
)
populations_list <- levels(population_labels$population)

# Functions
calculate_likelihood_ratio <- function(shared_alleles, genotype_match, pA, pB, k_values) {
  
  if (shared_alleles == 0) {
    return(k_values$k0)
  } else if (shared_alleles == 1) {
    
    Rxp <- switch(genotype_match,
                  "AA-AA" = pA,
                  "AA-AB" = 2 * pA,
                  "AB-AA" = 2 * pA,
                  "AB-AC" = 4 * pA,
                  "AB-AB" = (4 * pA * pB) / (pA + pB),
                  stop("Invalid genotype match for 1 shared allele.")
    )
    return(k_values$k0 + (k_values$k1 / Rxp))
    
  } else if (shared_alleles == 2) {
    Rxp <- switch(genotype_match,
                  "AA-AA" = pA,
                  "AB-AB" = (4 * pA * pB) / (pA + pB),
                  stop("Invalid genotype match for 2 shared alleles.")
    )
    Rxu <- ifelse(genotype_match == "AA-AA", pA^2, 2 * pA * pB)
    return(k_values$k0 + (k_values$k1 / Rxp) + (k_values$k2 / Rxu))
  } else {
    return(NA)
  }  
  
}

kinship_calculation <- function(row, allele_frequency_data, kinship_matrix) {
  alleles_ind1 <- as.character(c(row$ind1_allele1, row$ind1_allele2))
  alleles_ind2 <- as.character(c(row$ind2_allele1, row$ind2_allele2))
  
  # Get shared alleles and their counts
  shared_alleles_vector <- intersect(alleles_ind1, alleles_ind2)
  
  # Get unique alleles by excluding shared alleles
  unique_alleles_ind1 <- alleles_ind1[!alleles_ind1 %in% shared_alleles_vector]
  unique_alleles_ind2 <- alleles_ind2[!alleles_ind2 %in% shared_alleles_vector]
  
  # Create the allele_map using c() for combining alleles and setting names directly
  allele_map <- setNames(c(shared_alleles_vector, unique_alleles_ind1, unique_alleles_ind2), 
                         LETTERS[1:(length(shared_alleles_vector) + length(unique_alleles_ind1) + length(unique_alleles_ind2))])
  
  # Use match to label alleles more efficiently
  labeled_alleles_ind1 <- LETTERS[match(alleles_ind1, allele_map)]
  labeled_alleles_ind2 <- LETTERS[match(alleles_ind2, allele_map)]
  
  # Generate genotype strings
  genotype_ind1 <- paste(sort(labeled_alleles_ind1), collapse = "")
  genotype_ind2 <- paste(sort(labeled_alleles_ind2), collapse = "")
  genotype_match <- paste(genotype_ind1, genotype_ind2, sep = "-")
  
  # Count the shared alleles
  shared_alleles <- length(shared_alleles_vector)
  
  allele_freqs <- dplyr::filter(allele_frequency_data, population == row$population, marker == row$locus)
  pA <- allele_freqs$frequency[allele_freqs$allele == allele_map["A"]]
  pB <- allele_freqs$frequency[allele_freqs$allele == allele_map["B"]]
  
  kinship_calculations <- kinship_matrix[, .(relationship_known = row$relationship_type, 
                                             relationship_tested = relationship_type,
                                             LR = calculate_likelihood_ratio(shared_alleles, genotype_match, pA, pB, .SD))]
  
  cbind(row, kinship_calculations)
}

process_loci <- function(row, allele_frequency_data, kinship_matrix) {
  final_row <- kinship_calculation(row, allele_frequency_data, kinship_matrix)
  return(final_row)
}

process_individuals_genotypes <- function(individuals_genotypes, df_allelefreq, kinship_matrix) {
  # Map over rows of individuals_genotypes and process each row
  final_individuals_genotypes <- individuals_genotypes |>
    future_pmap(function(...) {
      row <- tibble(...)
      res <- process_loci(row, df_allelefreq, kinship_matrix)
      return(res)
    }, .progress = TRUE) |> 
       rbindlist()
  
  return(final_individuals_genotypes)
}

calculate_combined_lrs <- function(final_results, loci_lists) {
  final_results <- as.data.table(final_results)
  final_results[sapply(final_results, is.infinite)] <- NA
  combined_lrs <- final_results[, .(
    core_13 = prod(LR[locus %in% loci_lists$core_13], na.rm = TRUE),
    identifiler_15 = prod(LR[locus %in% loci_lists$identifiler_15], na.rm = TRUE),
    expanded_20 = prod(LR[locus %in% loci_lists$expanded_20], na.rm = TRUE),
    supplementary = prod(LR[locus %in% loci_lists$supplementary], na.rm = TRUE),
    autosomal_29 = prod(LR[locus %in% loci_lists$autosomal_29], na.rm = TRUE)
  ), by = .(population, relationship_known, relationship_tested, sim_id)]
  combined_lrs <- melt(combined_lrs,
                       id.vars = c("population", "relationship_known", "relationship_tested", "sim_id"),
                       measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                       variable.name = "loci_set", value.name = "LR")
  return(combined_lrs)
}

# Process the individuals genotypes
log_message("Processing individuals genotypes...")
processing_time <- system.time({
  processed_genotypes <- process_individuals_genotypes(individuals_genotypes, df_allelefreq, kinship_matrix)
})
log_message(paste("Processed individuals genotypes in", processing_time["elapsed"], "seconds."))

# Calculate combined likelihood ratios
log_message("Calculating combined likelihood ratios...")
combined_lrs_time <- system.time({
  combined_lrs <- calculate_combined_lrs(processed_genotypes, loci_lists)
})
log_message(paste("Calculated combined likelihood ratios in", combined_lrs_time["elapsed"], "seconds."))

# Save results to CSV
log_message("Saving results to CSV files...")
fwrite(processed_genotypes, output_file)
fwrite(combined_lrs, summary_output_file)

# Save timing log to CSV
timing_log_df <- as.data.frame(rbind(c("combined_lrs",combined_lrs_time), 
                                     c("genotype_processing",processing_time)))

fwrite(timing_log_df, timing_log_file)

log_message("Simulation processing completed.")
