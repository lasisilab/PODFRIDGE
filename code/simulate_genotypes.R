# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tibble)
  library(furrr)
  library(data.table)
  library(future)
  library(parallel)
  library(doParallel)
}))
args <- commandArgs(trailingOnly = TRUE)
print(args)
if(length(args)>4){
  use_remote_cluster<-as.numeric(args[5])
}
if(!exists("use_remote_cluster")){ #Add this parameter to the run script
  use_remote_cluster<-0
} #Use 0 if sending to external cluster (use 0 for tests on one machine)

# Set up cluster on one machine if required
if(use_remote_cluster==0){
  if( .Platform$OS.type == "windows" ){
    cl <- makeCluster(availableCores())  #specify how many cores to use
  } else { # use the fork cluster type on linux because its faster- not available for windows
    cl <- makeCluster(availableCores(),type="FORK")  #specify how many cores to use
  }
  # Ensure the cluster is stopped when the script exits
  # on.exit(parallel::stopCluster(cl))
  future::plan(future::cluster, workers = cl)
} else { #if sending to remote cluster(s)
  future::plan(future::cluster, workers = ("clustername1")) #use syntax "server.remote.org" in workers if using an online cluster
}

#plan() options:
#use 'multisession' to run in parallel in separate R sessions on the same machine
#use 'multicore' to run futures in parallel in forked processes on the same machine- Linux only
#use 'cluster' to run in parallel on one or more machines
#use 'sequential' to test in series (no parallel processing)

# up the limit of memory available to future per core
options('future.globals.maxSize' = 1014*1024^2)

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Helper function to log function timings
timing_log <- list()

log_function_time <- function(func, name, ...) {
  start_time <- Sys.time()
  result <- func(...)
  end_time <- Sys.time()
  duration <- as.numeric(difftime(end_time, start_time, units = "secs"))

  if (!name %in% names(timing_log)) {
    timing_log[[name]] <- list(total = 0, count = 0, min = Inf, max = -Inf, times = c())
  }

  timing_log[[name]]$total <- timing_log[[name]]$total + duration
  timing_log[[name]]$count <- timing_log[[name]]$count + 1
  timing_log[[name]]$min <- min(timing_log[[name]]$min, duration)
  timing_log[[name]]$max <- max(timing_log[[name]]$max, duration)
  timing_log[[name]]$times <- c(timing_log[[name]]$times, duration)

  return(list(result = result, timing_log = timing_log))
}

# Read Command-Line Arguments

n_sims_related <- as.numeric(args[1])
n_sims_unrelated <- as.numeric(args[2])
job_id <- as.character(args[3])

print("Creating directories")
print(getwd())
# Create output folder with job ID
output_dir <- file.path("data", "sims", paste0("simulation_", job_id))
dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "processed_genotypes.csv")
timing_log_file <- file.path(output_dir, "timing_log.csv")

# Log the start of the process
log_message("Starting genotype simulation...")

# Load Allele Frequencies Data
log_message("Loading allele frequencies data...")
allele_freq_time <- system.time({
  df_allelefreq <- fread("data/df_allelefreq_combined.csv")
  # df_allelefreq <- fread("data/syn_data.csv")
  df_allelefreq <- df_allelefreq[population != "all"] # Filter out "all" population
  df_allelefreq[, allele := as.character(allele)]
  df_allelefreq = as.data.table(df_allelefreq)

})
log_message(paste("Loaded allele frequencies data in", allele_freq_time["elapsed"], "seconds."))

# Extract unique loci
log_message("Extracting unique loci...")
loci_list <- unique(df_allelefreq$marker)

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

# Functions #generate_simulation_setup not run in parallel as already running quickly
generate_simulation_setup <- function(kinship_matrix, population_list, num_related, num_unrelated) {
  simulation_setup <- data.table(
    population = character(),
    relationship_type = character(),
    num_simulations = integer(),
    stringsAsFactors = FALSE
  )
  for (population in population_list) {
    for (relationship in kinship_matrix$relationship_type) {
      num_simulations <- ifelse(relationship == "unrelated", num_unrelated, num_related)
      l = list(simulation_setup, data.table(
        population = population,
        relationship_type = relationship,
        num_simulations = num_simulations
      ))
      simulation_setup <- rbindlist(l)
    }
  }
  return(simulation_setup)
}

initialize_individuals_pair <- function(population, relationship_type, sim_id, loci_list) {
  num_loci <- length(loci_list)
  nl <- character(num_loci)
  individuals_genotypes <- data.table(
    population = rep(population, num_loci),
    relationship_type = rep(relationship_type, num_loci),
    sim_id = rep(sim_id, num_loci),
    locus = loci_list,
    ind1_allele1 = nl,
    ind1_allele2 = nl,
    ind2_allele1 = nl,
    ind2_allele2 = nl)
  return(individuals_genotypes)
}

simulate_genotypes <- function(row, df_allelefreq, kinship_matrix) {
  population <- row$population
  locus <- row$locus
  relationship <- row$relationship_type

  allele_freqs <- df_allelefreq[which(df_allelefreq$population == population & df_allelefreq$marker == locus & df_allelefreq$frequency>0),]

  if (nrow(allele_freqs) == 0) {
    stop(paste("No valid alleles found for population ", population, " at locus", locus))
  }

  alleles <- allele_freqs$allele
  frequencies <- allele_freqs$frequency
  frequencies <- round(frequencies, 6)
 # print(summary(frequencies))

  ind1_alleles <- sample(alleles, size = 2, replace = TRUE, prob = frequencies)
  kinship_coeffs <- kinship_matrix[which(kinship_matrix$relationship_type == relationship), ]
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

process_individuals_genotypes <- function(sim_id,individuals_genotypes, df_allelefreq, kinship_matrix) {
  final_individuals_genotypes <- individuals_genotypes |>
    furrr::future_pmap(function(...) {
      res <- log_function_time(simulate_genotypes, "simulate_genotypes", list(...), df_allelefreq, kinship_matrix)
      return(res$result)
    }, seed = TRUE) |>
    data.table::rbindlist()
  final_individuals_genotypes$sim_id<- sim_id
  return(final_individuals_genotypes)
}

#     })
#   fwrite(final_results, gsub("processed","v1_processed",output_file))
#  })
#}

process_simulation_setup <- function(simulation_setup, df_allelefreq, kinship_matrix, loci_list, output_file) {

  process_time <- system.time({
    final_results <- simulation_setup |>
      future_pmap_dfr(function(population, relationship_type, num_simulations){
        cs<-seq(from=1,to=num_simulations,by=1)
        tibble(
          foreach(i = cs, .combine="rbind") %dopar% {
            sim_id = cs[i]
            individuals_genotypes <- initialize_individuals_pair(population, relationship_type, sim_id, loci_list)
            processed_genotypes <- log_function_time(process_individuals_genotypes, "process_individuals_genotypes", sim_id,individuals_genotypes, df_allelefreq, kinship_matrix)
            return(processed_genotypes$result)
          })
      }, .options = furrr::furrr_options(seed = NULL))

    fwrite(final_results, output_file,append=TRUE)
  })
  #  class(final_results)<-"data.table"

  log_message(paste("Processing completed in", process_time["elapsed"], "seconds."))

}

# Execute Simulation Setup and Processing (this is already too fast to benefit from being run in parallel)
setup_res <- log_function_time(generate_simulation_setup, "generate_simulation_setup", kinship_matrix, populations_list, n_sims_related, n_sims_unrelated)
simulation_setup <- setup_res$result
#simulation_setup is a table of all possible relationship and population combinations and the number of simulations to assign to each.
#Is this necessary as all unrelated will be assigned the same number of simulations etc

timing_log <- setup_res$timing_log
proc_res <- log_function_time(process_simulation_setup, "process_simulation_setup", simulation_setup, df_allelefreq, kinship_matrix, loci_list, output_file)
timing_log <- proc_res$timing_log

#proc_res_v1 <- log_function_time(process_simulation_setup_v1, "process_simulation_setup_v1", simulation_setup, df_allelefreq, kinship_matrix, loci_list, output_file)

#timing_log <- proc_res_v1$timing_log

log_message("Genotype simulation completed.")

# Save timing log to CSV
timing_log_df <- tibble(
  function_name = names(timing_log),
  total_time = sapply(timing_log, function(x) x$total),
  count = sapply(timing_log, function(x) x$count),
  min_time = sapply(timing_log, function(x) x$min),
  max_time = sapply(timing_log, function(x) x$max),
  avg_time = sapply(timing_log, function(x) x$total / x$count)
)

fwrite(timing_log_df, timing_log_file)
print(paste0("Timing log written to ",timing_log_file))
gc()
