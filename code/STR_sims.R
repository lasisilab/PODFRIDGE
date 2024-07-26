# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(furrr)
  library(progressr)
  library(data.table)
  library(ggplot2)
  library(future)
  library(parallel)
}))

# Set up cluster
cl <- makeCluster(availableCores())
plan(cluster, workers = cl)

# Ensure the cluster is stopped when the script exits
on.exit(parallel::stopCluster(cl))

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Helper function to log function timings
timing_log <- list()  # Initialize the timing_log list outside the function

log_function_time <- function(func, name, ...) {
  start_time <- Sys.time()
  result <- func(...)  # Call the wrapped function
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  # Update the timing_log
  if (!name %in% names(timing_log)) {
    timing_log[[name]] <<- list(total = 0, count = 0, min = Inf, max = -Inf, times = c())
  }
  timing_log[[name]]$total <<- timing_log[[name]]$total + duration
  timing_log[[name]]$count <<- timing_log[[name]]$count + 1
  timing_log[[name]]$min <<- min(timing_log[[name]]$min, duration)
  timing_log[[name]]$max <<- max(timing_log[[name]]$max, duration)
  timing_log[[name]]$times <<- c(timing_log[[name]]$times, duration)

  return(list(result = result, timings = timing_log)) # Return the result of func() and the timing_log
}

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
n_sims_related <- as.numeric(args[1])
n_sims_unrelated <- as.numeric(args[2])
output_dir <- args[3]  # Pass the output directory as an argument

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

# Define Kinship Matrix
kinship_matrix <- tibble(
  relationship_type = factor(
    c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"),
    levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated")
  ),
  k0 = c(0, 1/4, 1/2, 7/8, 15/16, 1),
  k1 = c(1, 1/2, 1/2, 1/8, 1/16, 0),
  k2 = c(0, 1/4, 0, 0, 0, 0)
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

generate_simulation_setup <- function(kinship_matrix, population_list, num_related, num_unrelated) {
  simulation_setup <- data.frame(
    population = character(),
    relationship_type = character(),
    num_simulations = integer(),
    stringsAsFactors = FALSE
  )
  for (population in population_list) {
    for (relationship in kinship_matrix$relationship_type) {
      num_simulations <- ifelse(relationship == "unrelated", num_unrelated, num_related)
      simulation_setup <- rbind(simulation_setup, data.frame(
        population = population,
        relationship_type = relationship,
        num_simulations = num_simulations
      ))
    }
  }
  return(simulation_setup)
}

initialize_individuals_pair <- function(population, relationship_type, sim_id, loci_list) {
  num_loci <- length(loci_list)
  individuals_genotypes <- data.table(
    population = rep(population, num_loci),
    relationship_type = rep(relationship_type, num_loci),
    sim_id = rep(sim_id, num_loci),
    locus = loci_list,
    ind1_allele1 = character(num_loci),
    ind1_allele2 = character(num_loci),
    ind2_allele1 = character(num_loci),
    ind2_allele2 = character(num_loci),
    shared_alleles = integer(num_loci),
    genotype_match = character(num_loci),
    LR = numeric(num_loci)
  )
  return(individuals_genotypes)
}

simulate_genotypes <- function(row, df_allelefreq, kinship_matrix) {
  population <- row$population
  locus <- row$locus
  relationship <- row$relationship_type

  allele_freqs <- df_allelefreq |>
    filter(population == !!population, marker == !!locus, frequency > 0)

  if (nrow(allele_freqs) == 0) {
    stop(paste("No valid alleles found for population", population, "and locus", locus))
  }

  alleles <- allele_freqs$allele
  frequencies <- allele_freqs$frequency

  frequencies <- round(frequencies / sum(frequencies), 6)
  valid_indices <- frequencies > 0
  alleles <- alleles[valid_indices]
  frequencies <- frequencies[valid_indices]

  ind1_alleles <- sample(alleles, size = 2, replace = TRUE, prob = frequencies)

  kinship_coeffs <- kinship_matrix[kinship_matrix$relationship_type == relationship, ]
  relationship_choice <- sample(c('none', 'one', 'both'), size = 1, prob = c(kinship_coeffs$k0, kinship_coeffs$k1, kinship_coeffs$k2))

  if (relationship_choice == 'none') {
    ind2_alleles <- sample(alleles, size = 2, replace = TRUE, prob = frequencies)
  } else if (relationship_choice == 'one') {
    shared_allele <- sample(ind1_alleles, size = 1)
    non_shared_allele <- sample(alleles, size = 1, replace = TRUE, prob = frequencies)
    if (runif(1) > 0.5) {
      ind2_alleles <- c(shared_allele, non_shared_allele)
    } else {
      ind2_alleles <- c(non_shared_allele, shared_allele)
    }
  } else if (relationship_choice == 'both') {
    ind2_alleles <- ind1_alleles
  }

  row$ind1_allele1 <- ind1_alleles[1]
  row$ind1_allele2 <- ind1_alleles[2]
  row$ind2_allele1 <- ind2_alleles[1]
  row$ind2_allele2 <- ind2_alleles[2]

  return(row)
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

  k_values <- kinship_matrix[kinship_matrix$relationship_type == row$relationship_type, ]
  LR <- calculate_likelihood_ratio(shared_alleles, genotype_match, pA, pB, k_values$k0, k_values$k1, k_values$k2)

  row$shared_alleles <- shared_alleles
  row$genotype_match <- genotype_match
  row$LR <- LR

  return(row)
}

process_loci <- function(row, allele_frequency_data, kinship_matrix) {
  simulated_row <- log_function_time(simulate_genotypes, "simulate_genotypes", row, allele_frequency_data, kinship_matrix)
  final_row <- log_function_time(kinship_calculation, "kinship_calculation", simulated_row, allele_frequency_data, kinship_matrix)
  return(final_row)
}

process_individuals_genotypes <- function(individuals_genotypes, df_allelefreq, kinship_matrix) {
  final_individuals_genotypes <- individuals_genotypes |>
    future_pmap(~ log_function_time(process_loci, "process_loci", list(...), df_allelefreq, kinship_matrix), seed = TRUE) |>
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
  ), by = .(population, relationship_type, sim_id)]
  combined_lrs <- melt(combined_lrs,
                       id.vars = c("population", "relationship_type", "sim_id"),
                       measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                       variable.name = "loci_set", value.name = "LR")
  return(combined_lrs)
}

plot_and_save_results <- function(combined_lrs) {
  # Ensure factor levels are set correctly for plotting
  combined_lrs$relationship_type <- factor(combined_lrs$relationship_type, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))

  ggplot(combined_lrs, aes(x = relationship_type, y = LR, fill = population, color = population)) +
    geom_boxplot() +
    facet_wrap(~ loci_set, scales = "fixed") +
    labs(
      title = "LR Distributions Across Populations and Relationship Types",
      x = "Relationship Type",
      y = "LR",
      fill = "Population",
      color = "Population"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_y_log10() +
    coord_flip()

  ggsave(file.path(output_dir, "sim_log_lr_panel_plot.png"), width = 12, height = 8)

  summary_stats <- combined_lrs |>
    group_by(relationship_type, population, loci_set) |>
    summarize(
      mean_LR = mean(LR),
      lower_95 = quantile(LR, 0.025),
      upper_95 = quantile(LR, 0.975)
    ) |>
    ungroup()

  ggplot(summary_stats, aes(x = loci_set, y = mean_LR, group = population, color = population)) +
    geom_line(size = 1) +
    facet_wrap(~ relationship_type, scales = "free_y", ncol = 2) +
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
    ) +
    scale_color_manual(values = c("AfAm" = "red", "Cauc" = "blue", "Hispanic" = "green", "Asian" = "purple"))

  ggsave(file.path(output_dir, "sim_line_chart_lr.png"), width = 14, height = 10)
}

plot_proportions_exceeding_cutoffs <- function(proportions_exceeding_cutoffs) {
  # Ensure factor levels are set correctly for plotting
  proportions_exceeding_cutoffs$relationship_type <- factor(proportions_exceeding_cutoffs$relationship_type, levels = c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"))
  proportions_exceeding_cutoffs$population <- factor(proportions_exceeding_cutoffs$population, levels = c("AfAm", "Cauc", "Hispanic", "Asian"))

  proportions_long <- proportions_exceeding_cutoffs |>
    pivot_longer(cols = starts_with("proportion_exceeding"),
                 names_to = "Cutoff_Type", values_to = "Proportion",
                 names_prefix = "proportion_exceeding_")

  proportions_long$Cutoff_Type <- factor(proportions_long$Cutoff_Type, levels = c("fixed", "1", "0_1", "0_01"),
                                         labels = c("Fixed Cutoff (1.00)", "1% FPR", "0.1% FPR", "0.01% FPR"))

  ggplot(proportions_long, aes(x = relationship_type, y = Proportion, fill = population, color = population)) +
    geom_bar(stat = "identity", position = position_dodge()) +
    facet_wrap(~ Cutoff_Type + loci_set, scales = "fixed") +
    labs(
      title = "Proportions Exceeding Likelihood Cut-offs Across Relationship Types and Loci Sets",
      x = "Relationship Type",
      y = "Proportion Exceeding Cut-off",
      fill = "Population",
      color = "Population"
    ) +
    theme_minimal() +
    theme(
      axis.text.x = element_text(angle = 45, hjust = 1)
    ) +
    scale_fill_manual(values = c("AfAm" = "red", "Cauc" = "blue", "Hispanic" = "green", "Asian" = "purple")) +
    coord_flip()

  ggsave(file.path(output_dir, "sim_proportions_exceeding_cutoffs_combined.png"), width = 12, height = 8)
}

# Function to calculate cut-off values for 1%, 0.1%, and 0.01% FPR
calculate_cutoffs <- function(input_df, fp_rates) {
  cutoffs <- input_df |>
    filter(relationship_type == "unrelated") |>
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
    group_by(population, relationship_type, loci_set) |>
    summarize(
      proportion_exceeding_fixed = sum(exceeds_fixed_cutoff, na.rm = TRUE) / n(),
      proportion_exceeding_1 = sum(exceeds_cutoff_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_1 = sum(exceeds_cutoff_0_1, na.rm = TRUE) / n(),
      proportion_exceeding_0_01 = sum(exceeds_cutoff_0_01, na.rm = TRUE) / n(),
      n_related = n(),
      .groups = 'drop'
    ) |>
    filter(relationship_type != "unrelated")
  return(proportions_exceeding)
}

process_simulation_setup <- function(simulation_setup, df_allelefreq, kinship_matrix, loci_list, loci_lists, output_file, summary_output_file) {
  log_message("Processing simulation setup...")
  process_time <- system.time({
    final_results <- simulation_setup |>
      future_pmap_dfr(function(population, relationship_type, num_simulations) {
        purrr::map_dfr(1:num_simulations, function(sim_id) {
          individuals_genotypes <- initialize_individuals_pair(population, relationship_type, sim_id, loci_list)
          processed_genotypes <- log_function_time(process_individuals_genotypes, "process_individuals_genotypes", individuals_genotypes, df_allelefreq, kinship_matrix)
          return(processed_genotypes)
        })
      }, .progress = TRUE)

    if ("seed" %in% colnames(final_results)) {
      final_results <- final_results |> select(-seed)
    }

    fwrite(final_results, output_file)
    combined_lrs <- log_function_time(calculate_combined_lrs, "calculate_combined_lrs", final_results, loci_lists)
    fwrite(combined_lrs, summary_output_file)

    # Create plots and save
    log_function_time(plot_and_save_results, "plot_and_save_results", combined_lrs)

    # Calculate and save cutoffs
    cutoffs <- log_function_time(calculate_cutoffs, "calculate_cutoffs", combined_lrs, c(1, 0.1, 0.01))
    fwrite(cutoffs, file.path(output_dir, "sim_cutoffs.csv"))

    proportions_exceeding_cutoffs <- log_function_time(calculate_proportions_exceeding_cutoffs, "calculate_proportions_exceeding_cutoffs", combined_lrs, cutoffs)
    fwrite(proportions_exceeding_cutoffs, file.path(output_dir, "sim_proportions_exceeding_cutoffs.csv"))

    # Convert population to factor for plotting
    proportions_exceeding_cutoffs$population <- factor(proportions_exceeding_cutoffs$population, levels = c("AfAm", "Cauc", "Hispanic", "Asian"))

    # Plot proportions exceeding cutoffs
    log_function_time(plot_proportions_exceeding_cutoffs, "plot_proportions_exceeding_cutoffs", proportions_exceeding_cutoffs)
  })
  log_message(paste("Processing completed in", process_time["elapsed"], "seconds."))
}

# Execute Simulation Setup and Processing
simulation_setup <- log_function_time(generate_simulation_setup, "generate_simulation_setup", kinship_matrix, populations_list, n_sims_related, n_sims_unrelated)
log_function_time(process_simulation_setup, "process_simulation_setup", simulation_setup, df_allelefreq, kinship_matrix, loci_list, loci_lists, output_file, summary_output_file)
log_message("Simulation setup and processing completed.")

# Save timing log to CSV
timing_log_df <- tibble(
  function_name = names(timing_log),
  total_time = sapply(timing_log, function(x) x$total),
  count = sapply(timing_log, function(x) x$count),
  min_time = sapply(timing_log, function(x) x$min),
  max_time = sapply(timing_log, function(x) x$max),
  avg_time = sapply(timing_log, function(x) x$total / x$count)
)

# Save timing log to CSV (modify this part based on the revised output)
timing_log_df <- bind_rows(lapply(timing_log, as.data.frame))
# Add the function names as a new column
timing_log_df$function_name <- names(timing_log)
write_csv(timing_log_df, timing_log_file)
