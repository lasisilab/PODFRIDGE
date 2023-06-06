# Load Required Libraries
library(tidyverse)
library(furrr)
library(progressr)

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
n_sims_unrelated <- as.numeric(args[1])
n_sims_related <- as.numeric(args[2])

# Load Allele Frequencies Data
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

# Define Functions
## Function: df-ibd-probs
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

## Function: simulate_STRpairs
simulate_STRpairs <- function(population, relationship_type, df_allelefreq, df_ibdprobs, n_sims=1) {

  markers <- unique(df_allelefreq$marker)

  allele_frequencies_by_marker <- df_allelefreq %>%
    filter(population == population) %>%
    split(.$marker)

  allele_frequencies_by_marker <- map(allele_frequencies_by_marker, ~.x %>% pull(frequency) %>% setNames(.x$allele))

  prob_shared_alleles <- df_ibdprobs %>%
    filter(relationship == relationship_type) %>%
    select(k0, k1, k2) %>%
    unlist() %>%
    as.numeric()
  non_zero_indices <- which(prob_shared_alleles != 0)

  k0 <- prob_shared_alleles[1]
  k1 <- prob_shared_alleles[2]
  k2 <- prob_shared_alleles[3]

  results <- lapply(markers, function(current_marker) {
    marker_results <- lapply(seq_len(n_sims), function(replicate_id) {
      individual1 <- setNames(
        lapply(markers, function(current_marker) {
          allele_frequencies <- allele_frequencies_by_marker[[current_marker]]
          return(sample(names(allele_frequencies), size = 2, replace = TRUE, prob = allele_frequencies))
        }), markers)

      individual2 <- setNames(
        lapply(markers, function(current_marker) {
          allele_frequencies <- allele_frequencies_by_marker[[current_marker]]
          num_shared_alleles <- sample(non_zero_indices - 1, size = 1, prob = prob_shared_alleles[non_zero_indices])
          alleles_from_individual1 <- sample(individual1[[current_marker]], size = num_shared_alleles)
          alleles_from_population <- sample(names(allele_frequencies), size = 2 - num_shared_alleles, replace = TRUE, prob = allele_frequencies)
          return(c(alleles_from_individual1, alleles_from_population))
        }), markers)

      ind1_alleles <- individual1[[current_marker]]
      ind2_alleles <- individual2[[current_marker]]
      shared_alleles <- intersect(ind1_alleles, ind2_alleles)
      num_shared_alleles <- length(shared_alleles)

      R_Xp <- sum(purrr::map_dbl(shared_alleles, function(x) unlist(allele_frequencies_by_marker[[current_marker]][x])))
      R_Xu <- sum(purrr::map_dbl(ind1_alleles, function(x) unlist(allele_frequencies_by_marker[[current_marker]][x])) * purrr::map_dbl(ind2_alleles, function(x) unlist(allele_frequencies_by_marker[[current_marker]][x])))


      R <- k0
      if (R_Xp != 0) { R <- R + (k1 / R_Xp) }
      if (R_Xu != 0) { R <- R + (k2 / R_Xu) }

      log_R <- log(R)

      # Add the replicate_id column to the output tibble
      return(tibble(population = population,
                    relationship_type = relationship_type,
                    marker = current_marker,
                    num_shared_alleles = num_shared_alleles,
                    log_R = log_R,
                    replicate_id = replicate_id)) # Add this line
    })

    marker_results <- bind_rows(marker_results)
    return(marker_results)
  })

  result <- bind_rows(results)

  # Aggregate results for each replicate, summing num_shared_alleles and log_R values.
  result_by_replicate <- result %>%
    group_by(population, relationship_type, replicate_id) %>%
    summarise(num_shared_alleles_sum = sum(num_shared_alleles),
              log_R_sum = sum(log_R),
              .groups = "drop")

  return(result_by_replicate)

}


## Function: simulation_combinations
simulation_combinations <- function(df, n_sims_unrelated, n_sims_related) {
  # Define the list of relationship types
  relationship_types <- c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated')

  # Get unique populations from the input dataframe
  unique_populations <- unique(df$population)
  filtered_populations <- unique_populations[unique_populations != "all"]

  # Create a dataframe of all combinations of populations and relationship types
  combinations <- expand_grid(population = filtered_populations, relationship_type = relationship_types)

  # Add the number of simulations for unrelated or related relationships
  combinations$n_sims <- ifelse(combinations$relationship_type == "unrelated", n_sims_unrelated, n_sims_related)

  return(combinations)
}

# Function to capitalize the first letter of a string
ucfirst <- function(s) {
  paste(toupper(substring(s, 1,1)), substring(s, 2), sep = "")
}

create_plot <- function(df, variable_to_plot, relationship_col, population_col) {
  # Set the population_shape variable as factor
  df$population_shape <- factor(df[[population_col]])

  # Create the plot
  p <- ggplot(df, aes(x = .data[[relationship_col]], y = .data[[variable_to_plot]],
                      color = .data[[population_col]], shape = .data[[population_col]], fill = .data[[population_col]])) +
    geom_boxplot(alpha = 0.4, position = position_dodge(width = 0.75)) +
    geom_point(position = position_jitterdodge(jitter.width = 0.1, dodge.width = 0.75), size = 1, alpha = 0.6) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1)) +
    scale_x_discrete(limits = c('parent_child', 'full_siblings', 'half_siblings', 'cousins', 'second_cousins', 'unrelated')) +
    scale_shape_manual(values = c(16, 17, 15, 18)) + # Change these values to desired shapes
    labs(title = paste(ucfirst(variable_to_plot), "by Relationship Type and Population"),
         x = "Relationship Type",
         y = ucfirst(variable_to_plot),
         color = "Population",
         shape = "Population",
         fill = "Population")

  # Save the plot to the /output folder with a custom file name
  save_plot <- function(plot, plot_name) {
    ggsave(filename = paste0("output/", plot_name, ".png"), plot = plot, height = 6, width = 8, units = "in")
  }

  # Call the save_plot function to save the plot
  save_plot(p, paste("plot_", variable_to_plot, sep = ""))

  return(p)
}


# Execute main script pipeline
result_combinations <- simulation_combinations(df_freq, n_sims_unrelated = n_sims_unrelated, n_sims_related = n_sims_related)

plan(multisession)
handlers(global = TRUE)
handlers("progress")

# Set future.rng.onMisuse to "ignore" to suppress the warnings
options(future.rng.onMisuse = "ignore")

results_parallel <- with_progress({
  result_combinations %>%
    future_pmap_dfr(function(population, relationship_type, n_sims) {
      future.seed <- sample.int(.Machine$integer.max, 1) # Randomly assign a seed value at each iteration
      simulate_STRpairs(population, relationship_type, df_allelefreq = df_freq, df_ibdprobs = df_ibdprobs, n_sims = n_sims)
    }, .progress = TRUE)
})

# Save result as CSV in the /data folder
write_csv(results_parallel, file = "data/simulation_results.csv")

# Create Plots
df_plt_final <- results_parallel %>% select(-replicate_id)
p_log_R_sum <- create_plot(df_plt_final, "log_R_sum", "relationship_type", "population")
p_num_shared_alleles_sum <- create_plot(df_plt_final, "num_shared_alleles_sum", "relationship_type", "population")
