library(tidyverse)
library(furrr)
library(progressr)

args <- commandArgs(trailingOnly = TRUE)
n_sims_unrelated <- as.numeric(args[1])
n_sims_related <- as.numeric(args[2])

marker_list <- readLines("data/marker_list.txt")

file_paths <- list.files(path = "data", pattern = "1036_.*\\.csv", full.names = TRUE)
df_list <- lapply(file_paths, function(path) {
  read_csv(path, col_types = cols(
    marker = col_character(),
    allele = col_double(),
    frequency = col_double(),
    population = col_character()
  ))
})
df_freq <- bind_rows(df_list)

# Filter the data frame based on the marker list
df_freq <- df_freq %>% filter(marker %in% marker_list)

create_df_ibdprobs <- function() {
  tibble(
    relationship =
      c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins", "unrelated"),
    k0 = c(0, 1/4, 1/2, 7/8, 15/16, 1),
    k1 = c(1, 1/2, 1/2, 1/8, 1/16, 0),
    k2 = c(0, 1/4, 0, 0, 0, 0)
  )
}
df_ibdprobs <- create_df_ibdprobs()

get_allele_frequencies <- function(population, marker, num_shared_alleles, A, B, df_allelefreq) {
  # Initialize the output values for p_A and p_B
  p_A <- NA
  p_B <- NA

  # Check whether there are shared alleles
  if (num_shared_alleles > 0) {
    # Get the frequency of allele A in the df_allelefreq dataset
    p_A <- df_allelefreq %>%
      filter(population == population, marker == marker, allele == A) %>%
      group_by(population, marker, allele) %>%
      pull(frequency) %>%
      first()

    # If two alleles are shared and B is not NA, get the frequency of allele B
    if (num_shared_alleles == 2 && !is.na(B)) {
      p_B <- df_allelefreq %>%
        filter(population == population, marker == marker, allele == B) %>%
        group_by(population, marker, allele) %>%
        pull(frequency) %>%
        first()
    }
  }

  # Return named list of p_A and p_B
  return(list(p_A = p_A, p_B = p_B))
}


# Define a function to calculate R_Xp values based on ind1_geno (c), ind2_geno (Q), p_A, and p_B values
get_R_Xp <- function(c, Q, p_A, p_B) {
  if (c == "AA") {
    if (Q == "AA") {
      return(p_A)
    } else if (Q == "AB") {
      return(2 * p_A)
    }
  } else if (c == "AB") {
    if (Q == "AA") {
      return(2 * p_A)
    } else if (Q == "AB") {
      return(4 * p_A * p_B / (p_A + p_B))
    } else if (Q == "AC") {
      return(4 * p_A)
    }
  }
  return(0)
}

# Define a function to calculate R_Xu values based on ind1_geno (c), ind2_geno (Q), p_A, and p_B values
get_R_Xu <- function(c, Q, p_A, p_B) {
  if (c == Q) {
    if (c == "AA") {
      return(p_A^2)
    } else if (c == "AB") {
      return(2 * p_A * p_B)
    }
  }
  return(0)
}

# Define a function to calculate R values using the relationship, R_Xp, and R_Xu values, and the global df_ibdprobs dataframe
calc_R <- function(relationship_type, R_Xp, R_Xu, df_ibdprobs) {

  # Get the corresponding kappa values based on the given relationship type
  df_ibdprobs$relationship <- as.character(df_ibdprobs$relationship)

  k_values <- df_ibdprobs %>%
    filter(relationship == relationship_type) %>%
    select(k0, k1, k2) %>%
    as.list()

  k0 <- k_values$k0
  k1 <- k_values$k1
  k2 <- k_values$k2

  # Calculate the R value using the provided formula with correct parenthesis and handling of R_Xp and R_Xu equal to zero
  R <- k0
  if (R_Xp != 0) {
    R <- R + (k1 / R_Xp)
  }
  if (R_Xu != 0) {
    R <- R + (k2 / R_Xu)
  }

  # Calculate the log(R) value
  log_R <- log(R)

  # Return R and log(R) as a list
  return(list(R = R, log_R = log_R))
}


shared_alleles <- function(ind_1_allele_1, ind_1_allele_2, ind_2_allele_1, ind_2_allele_2, population, marker, relationship_type, df_allelefreq, df_ibdprobs) {

  # Create allele vectors for individuals 1 and 2
  alleles_ind1 <- c(ind_1_allele_1, ind_1_allele_2)
  alleles_ind2 <- c(ind_2_allele_1, ind_2_allele_2)

  # Initialize variables to keep track of shared alleles and used alleles
  shared <- character()
  used_ind1 <- rep(FALSE, 2)
  used_ind2 <- rep(FALSE, 2)

  # Iterate through alleles of individual 1
  for (i in 1:2) {
    allele_ind1 <- alleles_ind1[i]

    # Iterate through alleles of individual 2
    for (j in 1:2) {
      # Check if alleles are not used and are the same
      if (!used_ind2[j] && allele_ind1 == alleles_ind2[j]) {
        # Add the shared allele to the shared vector
        shared <- c(shared, allele_ind1)
        # Mark the alleles as used
        used_ind1[i] <- TRUE
        used_ind2[j] <- TRUE
        break
      }
    }
  }

  # Calculate the number of shared alleles and find unique shared alleles
  num_shared_alleles <- length(shared)
  unique_shared <- unique(shared)

  # Create a letter map for shared and non-shared alleles using capitalized letters
  letter_map <- setNames(LETTERS[seq_along(unique_shared)], unique_shared)
  remaining_alleles <- union(alleles_ind1, alleles_ind2)
  non_shared <- setdiff(remaining_alleles, unique_shared)
  letter_map <- c(letter_map, setNames(LETTERS[length(unique_shared) + seq_along(non_shared)], non_shared))

  # Assign capitalized letters to alleles based on the letter map and create sorted vectors
  ind1_geno <- sort(letter_map[alleles_ind1])
  ind2_geno <- sort(letter_map[alleles_ind2])

  # Initialize A and B
  A <- NA_character_
  B <- NA_character_

  if (num_shared_alleles > 0) {
    A <- names(letter_map)[1]

    if (length(unique_shared) > 1) {
      B <- names(letter_map)[2]
    }
  }

  # Sort ind1_geno and ind2_geno alphabetically and paste them together as a string vector
  ind1_geno <- paste0(sort(ind1_geno), collapse = "")
  ind2_geno <- paste0(sort(ind2_geno), collapse = "")

  # Return a single-row tibble directly
  result <- tibble::tibble(
    num_shared_alleles = num_shared_alleles,
    ind1_geno = ind1_geno,
    ind2_geno = ind2_geno,
    A = A,
    B = B
  )

  # Call the get_allele_frequencies function with required arguments
  allele_freqs <- get_allele_frequencies(population, marker, num_shared_alleles, A, B, df_allelefreq)

  # Add p_A and p_B values to the result tibble
  result$p_A <- allele_freqs$p_A
  result$p_B <- allele_freqs$p_B

  # Call the new functions to calculate R_Xp and R_Xu values and add them to the result tibble
  result$R_Xp <- get_R_Xp(result$ind1_geno, result$ind2_geno, result$p_A, result$p_B)
  result$R_Xu <- get_R_Xu(result$ind1_geno, result$ind2_geno, result$p_A, result$p_B)

  # Call the new calc_R function to calculate R and log(R) values and add them to the result tibble
  R_values <- calc_R(relationship_type, result$R_Xp, result$R_Xu, df_ibdprobs)
  result$R <- R_values$R
  result$log_R <- R_values$log_R

  return(result)
}

# Define a function called 'simulate_STRpairs' to generate pairs of STR profiles for a known relationship type and calculate relatedness for tested relationship type
simulate_STRpairs <- function(population, known_relationship_type, tested_relationship_type, df_allelefreq, df_ibdprobs, n_sims=1) {

  # Extract unique markers from the df_allelefreq dataframe
  markers <- unique(df_allelefreq$marker)

  # Filter df_allelefreq by population, then split by marker and extract frequencies
  allele_frequencies_by_marker <- df_allelefreq %>%
    filter(population == population) %>%
    split(.$marker)

  allele_frequencies_by_marker <- map(allele_frequencies_by_marker, ~.x %>% pull(frequency) %>% setNames(.x$allele))

  # Get probabilities of shared alleles (k0, k1, k2) for the known_relationship_type from the df_ibdprobs dataframe
  prob_shared_alleles <- df_ibdprobs %>%
    # Filter df_ibdprobs to keep rows with relationship equal to known_relationship_type
    filter(relationship == known_relationship_type) %>%
    # Select k0, k1, and k2 columns (probabilities of sharing 0, 1, or 2 alleles identical by descent)
    select(k0, k1, k2) %>%
    # Convert the dataframe columns into a single numeric vector
    unlist() %>%
    # Ensure the output vector has a numeric datatype
    as.numeric()

  # Identify the indices of non-zero probabilities in the prob_shared_alleles vector
  non_zero_indices_known <- which(prob_shared_alleles != 0)


  # Assign k0, k1, and k2 probabilities for known_relationship_type
  k0_known <- prob_shared_alleles[1]
  k1_known <- prob_shared_alleles[2]
  k2_known <- prob_shared_alleles[3]

  # Get probabilities of shared alleles (k0, k1, k2) for the tested_relationship_type from the df_ibdprobs dataframe
  prob_shared_alleles_tested <- df_ibdprobs %>%
    # Filter df_ibdprobs to keep rows with relationship equal to tested_relationship_type
    filter(relationship == tested_relationship_type) %>%
    # Select k0, k1, and k2 columns (probabilities of sharing 0, 1, or 2 alleles identical by descent)
    select(k0, k1, k2) %>%
    # Convert the dataframe columns into a single numeric vector
    unlist() %>%
    # Ensure the output vector has a numeric datatype
    as.numeric()

  # Identify the indices of non-zero probabilities in the prob_shared_alleles_tested vector
  non_zero_indices_tested <- which(prob_shared_alleles_tested != 0)


  # Assign k0, k1, and k2 probabilities for tested_relationship_type
  k0_tested <- prob_shared_alleles_tested[1]
  k1_tested <- prob_shared_alleles_tested[2]
  k2_tested <- prob_shared_alleles_tested[3]

  # Run simulations for each marker
  results <- lapply(markers, function(current_marker) {
    marker_results <- lapply(seq_len(n_sims), function(replicate_id) {

      # Generate individual1's alleles based on the alleles frequencies
      individual1 <- setNames(
        # Iterate through the markers using the lapply function
        lapply(markers, function(current_marker) {
          # Retrieve the allele frequencies for the current_marker from the allele_frequencies_by_marker list
          allele_frequencies <- allele_frequencies_by_marker[[current_marker]]
          # Sample two alleles according to their frequencies, replace = TRUE allows sampling the same allele twice
          return(sample(names(allele_frequencies), size = 2, replace = TRUE, prob = allele_frequencies))
          # Set the names of the list elements as the marker names
        }), markers)


      # Generate individual2's alleles based on the alleles frequencies and the number of shared alleles with individual1
      individual2 <- setNames(
        lapply(markers, function(current_marker) {
          # Retrieve the allele frequencies for the current_marker from the allele_frequencies_by_marker list
          allele_frequencies <- allele_frequencies_by_marker[[current_marker]]
          # Sample the number of shared alleles following the known relationship type probabilities, using non_zero_indices_known
          num_shared_alleles <- sample(non_zero_indices_known - 1, size = 1, prob = prob_shared_alleles[non_zero_indices_known])
          # Select the shared alleles from individual1
          alleles_from_individual1 <- sample(individual1[[current_marker]], size = num_shared_alleles)
          # Sample the remaining alleles from the population for individual2
          alleles_from_population <- sample(names(allele_frequencies), size = 2 - num_shared_alleles, replace = TRUE, prob = allele_frequencies)
          # Combine the shared alleles from individual1 and the remaining alleles from thepopulation
          return(c(alleles_from_individual1, alleles_from_population))
          # Set the names of the list elements as the marker names
        }), markers)

      # Calculate the number of shared alleles between the two individuals for the current marker
      ind1_alleles <- individual1[[current_marker]]
      ind2_alleles <- individual2[[current_marker]]

      shared_allele_result <- shared_alleles(ind1_alleles[1], ind1_alleles[2], ind2_alleles[1], ind2_alleles[2], population, current_marker, tested_relationship_type, df_allelefreq, df_ibdprobs)

      # Extract num_shared_alleles and log_R from the shared_allele_result:
      num_shared_alleles <- shared_allele_result$num_shared_alleles
      log_R <- shared_allele_result$log_R


      # Store the results from the current simulation in a tibble
      return(tibble(population = population,
                    known_relationship_type = known_relationship_type,
                    tested_relationship_type = tested_relationship_type,
                    k0_tested = k0_tested,
                    k1_tested = k1_tested,
                    k2_tested = k2_tested,
                    marker = current_marker,
                    num_shared_alleles = num_shared_alleles,
                    log_R = log_R,
                    replicate_id = replicate_id))
    })

    # Combine the results from the current marker
    marker_results <- bind_rows(marker_results)
    return(marker_results)
  })

  # Combine the results from all markers
  result <- bind_rows(results)

  # Summarize the results by replicate
  result_by_replicate <- result %>%
    group_by(population, known_relationship_type, tested_relationship_type, replicate_id) %>%
    summarise(num_shared_alleles_sum = sum(num_shared_alleles),
              log_R_sum = sum(log_R),
              .groups = "drop")

  return(result_by_replicate)

}


simulation_combinations_all_relationships <- function(df, n_sims_unrelated, n_sims_related) {
  unique_populations <- unique(df$population)
  filtered_populations <- unique_populations[unique_populations != "all"]

  relationship_types <- c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated')

  combinations <- expand_grid(population = filtered_populations, known_relationship_type = relationship_types, tested_relationship_type = relationship_types)

  combinations$n_sims <- ifelse(combinations$known_relationship_type == "unrelated", n_sims_unrelated, n_sims_related)

  return(combinations)
}

ucfirst <- function(s) {
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep = "")
}

create_plot <- function(df, variable_to_plot, known_relationship_col, tested_relationship_col, population_col) {
  df$population_shape <- factor(df[[population_col]])
  df[[known_relationship_col]] <- factor(df[[known_relationship_col]], levels = c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated'))

  p <- ggplot(df, aes(x = .data[[tested_relationship_col]], y = .data[[variable_to_plot]],
                      color = .data[[population_col]], shape = .data[[population_col]], fill = .data[[population_col]])) +
    geom_boxplot(alpha = 0.4, position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 1, alpha = 0.6) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated')) +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    labs(title = paste(ucfirst(variable_to_plot), "by Tested Relationship Type and Population"),
         x = "Tested Relationship Type",
         y = ucfirst(variable_to_plot),
         color = "Population",
         shape = "Population",
         fill = "Population") +
    facet_wrap(~ .data[[known_relationship_col]], ncol = 1)

  return(p)
}


save_plot <- function(plot, plot_name) {
  ggsave(filename = paste0("output/known_vs_tested_", plot_name, ".png"), plot = plot, height = 15, width = 8, units = "in")
}

create_heatmap <- function(df, known_relationship_col, tested_relationship_col) {
  df_mean_diff <- df %>%
    group_by(.data[[known_relationship_col]], .data[[tested_relationship_col]]) %>%
    summarize(mean_diff = mean(log_R_sum), .groups = "drop")

  relationship_order <- c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated')

  heatmap_plot <- ggplot(data = df_mean_diff, aes(x = .data[[tested_relationship_col]], y = .data[[known_relationship_col]], fill = mean_diff)) +
    geom_tile(color = "black") +
    scale_fill_gradient2(low = "blue", mid = "white", high = "red", midpoint = 0) +
    scale_x_discrete(limits = relationship_order) +
    scale_y_discrete(limits = relationship_order) +
    theme_minimal() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    labs(title = "Heatmap of Mean Difference in Log_R_sum",
         x = "Tested Relationship Type",
         y = "Known Relationship Type",
         fill = "Mean Difference")

  return(heatmap_plot)
}



result_combinations <- simulation_combinations_all_relationships(df_freq, n_sims_unrelated = n_sims_unrelated, n_sims_related = n_sims_related)

plan(multisession, workers = availableCores()-1)
handlers(global = TRUE)
handlers("progress")

options(future.rng.onMisuse = "ignore")

results_parallel <- with_progress({
  result_combinations %>%
    future_pmap_dfr(function(population, known_relationship_type, tested_relationship_type, n_sims) {
      future.seed <- sample.int(.Machine$integer.max, 1)
      simulate_STRpairs(population, known_relationship_type, tested_relationship_type, df_allelefreq = df_freq, df_ibdprobs = df_ibdprobs, n_sims = n_sims)
    }, .progress = TRUE)
})

write_csv(results_parallel, file = "data/known_vs_tested_simulation_results.csv")

df_plt_final <- results_parallel %>% select(-replicate_id)
p_log_R_sum <- create_plot(df_plt_final, "log_R_sum", "known_relationship_type", "tested_relationship_type", "population")
p_num_shared_alleles_sum <- create_plot(df_plt_final, "num_shared_alleles_sum", "known_relationship_type", "tested_relationship_type", "population")

save_plot(p_log_R_sum, "plot_log_R_sum")
save_plot(p_num_shared_alleles_sum, "plot_num_shared_alleles_sum")

heatmap_log_R_sum <- create_heatmap(df_plt_final, "known_relationship_type", "tested_relationship_type")
ggsave(filename = "output/known_vs_tested_heatmap_log_R_sum.png", plot = heatmap_log_R_sum, height = 6, width = 8, units = "in")

create_plot_v2 <- function(df, variable_to_plot, known_relationship_col, tested_relationship_col, population_col) {
  df$population_shape <- factor(df[[population_col]])
  df[[tested_relationship_col]] <- factor(df[[tested_relationship_col]], levels = c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated'))

  p <- ggplot(df, aes(x = .data[[known_relationship_col]], y = .data[[variable_to_plot]],
                      color = .data[[population_col]], shape = .data[[population_col]], fill = .data[[population_col]])) +
    geom_boxplot(alpha = 0.4, position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 1, alpha = 0.6) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated')) +
    scale_shape_manual(values = c(16, 17, 15, 18)) +
    labs(title = paste(ucfirst(variable_to_plot), "by Tested Relationship Type and Population"),
         x = "Known Relationship Type",
         y = ucfirst(variable_to_plot),
         color = "Population",
         shape = "Population",
         fill = "Population") +
    facet_wrap(~ .data[[tested_relationship_col]], ncol = 1)

  return(p)
}

df_plt_final_v2 <- results_parallel %>% select(-replicate_id)

p_log_R_sum_v2 <- create_plot_v2(df_plt_final_v2, "log_R_sum", "known_relationship_type", "tested_relationship_type", "population")
p_num_shared_alleles_sum_v2 <- create_plot_v2(df_plt_final_v2, "num_shared_alleles_sum", "known_relationship_type", "tested_relationship_type", "population")

save_plot(p_log_R_sum_v2, "plot_log_R_sum_v2")
save_plot(p_num_shared_alleles_sum_v2, "plot_num_shared_alleles_sum_v2")




