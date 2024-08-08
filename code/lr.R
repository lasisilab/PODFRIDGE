# Load Required Libraries
library(tidyverse)
library(furrr)
library(data.table)
library(future)
library(parallel)

# Set up cluster
cl <- makeCluster(availableCores())
plan(cluster, workers = cl)

# Ensure the cluster is stopped when the script exits
on.exit(parallel::stopCluster(cl))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
slurm_job_id <- as.character(args[1])

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
loci_list <- df_allelefreq |> 
  pull(marker) |> 
  unique()

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
individuals_genotypes <- read.csv("data/sims/processed_genotypes.csv")

# Define Kinship Matrix
kinship_matrix <- tibble(
  relationship_type = factor(
    c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins"),
    levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins")
  ),
  k0 = c(0, 1/4, 1/2, 7/8, 15/16),
  k1 = c(1, 1/2, 1/2, 1/8, 1/16),
  k2 = c(0, 1/4, 0, 0, 0)
)

# Define Populations
population_labels <- tibble(
  population = factor(
    c("AfAm", "Cauc", "Hispanic", "Asian"),
    levels = c("AfAm", "Cauc", "Hispanic", "Asian")
  ),
  label = c("African American", "Caucasian", "Hispanic", "Asian")
)
populations_list <- levels(population_labels$population)

# Functions
calculate_likelihood_ratio <- function(shared_alleles, genotype_match = NULL, pA = NULL, pB = NULL, k0, k1, k2) {
  if (shared_alleles == 0) {
    LR <- k0
    return(LR)
  }
  if (shared_alleles == 1) {
    if (genotype_match == "AA-AA") {
      Rxp <- pA
    } else if (genotype_match == "AA-AB" | genotype_match == "AB-AA") {
      Rxp <- 2 * pA
    } else if (genotype_match == "AB-AC") {
      Rxp <- 4 * pA
    } else if (genotype_match == "AB-AB") {
      Rxp <- (4 * (pA * pB)) / (pA + pB)
    } else {
      stop("Invalid genotype match for 1 shared allele.")
    }
    LR <- k0 + (k1 / Rxp)
    return(LR)
  }
  if (shared_alleles == 2) {
    if (genotype_match == "AA-AA") {
      Rxp <- pA
      Rxu <- pA^2
    } else if (genotype_match == "AB-AB") {
      Rxp <- (4 * pA * pB) / (pA + pB)
      Rxu <- 2 * pA * pB
    } else {
      stop("Invalid genotype match for 2 shared alleles.")
    }
    LR <- k0 + (k1 / Rxp) + (k2 / Rxu)
    return(LR)
  }
}

kinship_calculation <- function(row, allele_frequency_data, kinship_matrix) {
  alleles_ind1 <- as.character(c(row$ind1_allele1, row$ind1_allele2))
  alleles_ind2 <- as.character(c(row$ind2_allele1, row$ind2_allele2))
  
  shared_alleles_vector <- intersect(alleles_ind1, alleles_ind2)
  unique_alleles_ind1 <- setdiff(alleles_ind1, shared_alleles_vector)
  unique_alleles_ind2 <- setdiff(alleles_ind2, shared_alleles_vector)
  
  allele_map <- list()
  next_label <- 1
  
  for (allele in shared_alleles_vector) {
    allele_map[[LETTERS[next_label]]] <- allele
    next_label <- next_label + 1
  }
  
  for (allele in unique_alleles_ind1) {
    if (!(allele %in% allele_map)) {
      allele_map[[LETTERS[next_label]]] <- allele
      next_label <- next_label + 1
    }
  }
  
  for (allele in unique_alleles_ind2) {
    if (!(allele %in% allele_map)) {
      allele_map[[LETTERS[next_label]]] <- allele
      next_label <- next_label + 1
    }
  }
  
  allele_map <- unlist(allele_map)
  labeled_alleles_ind1 <- sapply(as.character(alleles_ind1), function(x) names(allele_map)[which(allele_map == x)])
  labeled_alleles_ind2 <- sapply(as.character(alleles_ind2), function(x) names(allele_map)[which(allele_map == x)])
  
  shared_alleles <- length(shared_alleles_vector)
  genotype_ind1 <- paste(sort(labeled_alleles_ind1), collapse = "")
  genotype_ind2 <- paste(sort(labeled_alleles_ind2), collapse = "")
  genotype_match <- paste(genotype_ind1, genotype_ind2, sep = "-")
  
  allele_freqs <- dplyr::filter(allele_frequency_data, population == row$population, marker == row$locus)
  if (nrow(allele_freqs) == 0) {
    stop("No allele frequencies found for the given population and locus.")
  }
  
  A_allele <- ifelse("A" %in% names(allele_map), allele_map[["A"]], NA)
  B_allele <- ifelse("B" %in% names(allele_map), allele_map[["B"]], NA)
  
  pA <- ifelse(any(allele_freqs$allele == A_allele), allele_freqs$frequency[allele_freqs$allele == A_allele], NA)
  pB <- ifelse(any(allele_freqs$allele == B_allele), allele_freqs$frequency[allele_freqs$allele == B_allele], NA)
  
  if (is.na(pA)) {
    stop("Allele frequency for A is missing.")
  }
  if (is.na(pB) && length(shared_alleles_vector) > 1) {
    stop("Allele frequency for B is missing.")
  }
  
  kinship_calculations <- lapply(kinship_matrix$relationship_type, function(rel_type) {
    k_values <- kinship_matrix[kinship_matrix$relationship_type == rel_type, ]
    LR <- calculate_likelihood_ratio(shared_alleles, genotype_match, pA, pB, k_values$k0, k_values$k1, k_values$k2)
    data.table(
      relationship_known = row$relationship_type,
      relationship_tested = rel_type,
      LR = LR
    )
  })
  
  kinship_calculations <- rbindlist(kinship_calculations)
  row <- cbind(row, kinship_calculations)
  return(row)
}

process_loci <- function(row, allele_frequency_data, kinship_matrix) {
  final_row <- kinship_calculation(row, allele_frequency_data, kinship_matrix)
  return(final_row)
}

process_individuals_genotypes <- function(individuals_genotypes, df_allelefreq, kinship_matrix) {
  # Map over rows of individuals_genotypes and process each row
  final_individuals_genotypes <- individuals_genotypes %>%
    future_pmap(function(...) {
      row <- tibble(...)
      res <- process_loci(row, df_allelefreq, kinship_matrix)
      return(res)
    }, .progress = TRUE) %>% 
    bind_rows()
  
  return(final_individuals_genotypes)
}

calculate_combined_lrs <- function(final_results, loci_lists) {
  final_results <- as.data.table(final_results)
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
write.csv(processed_genotypes, output_file, row.names = FALSE)
write.csv(combined_lrs, summary_output_file, row.names = FALSE)

# Save timing log to CSV
timing_log_df <- as.data.frame(rbind(c("combined_lrs",combined_lrs_time), 
                                     c("genotype_processing",processing_time)))

write_csv(timing_log_df, timing_log_file)

log_message("Simulation processing completed.")
