# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(furrr)
  library(data.table)
  library(future)
  library(parallel)
  library(doParallel)
  library(rslurm)
}))

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
n_sims_related <- as.numeric(args[1])
n_sims_unrelated <- as.numeric(args[2])
job_id <- as.character(args[3])

if(length(args)>3){
  use_remote_cluster<-as.numeric(args[4])
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
  on.exit(parallel::stopCluster(cl))
  plan(cluster, workers = cl)
} else { #if sending to remote cluster(s)
  plan(cluster, workers = ("clustername1")) #use syntax "server.remote.org" in workers if using an online cluster
}

#plan() options:
#use 'multisession' to run in parallel in separate R sessions on the same machine
#use 'multicore' to run futures in parallel in forked processes on the same machine- Linux only
#use 'cluster' to run in parallel on one or more machines
#use 'sequential' to test in series (no parallel processing)

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
#  df_allelefreq <- fread("data/syn_data.csv")
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

  allele_freqs <- df_allelefreq[which(df_allelefreq$population == population & df_allelefreq$marker == locus),]

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

#purrr::map inside the future_pmap_dfr loop seems to slow this section (solution below is faster), but this version also works.
#process_simulation_setup_v1 <- function(simulation_setup, df_allelefreq, kinship_matrix, loci_list, output_file) {
 #   process_time <- system.time({
  #  class(simulation_setup)<-"data.frame"
   # final_results <- simulation_setup |>
    # furrr::future_pmap_dfr(.f=function(population, relationship_type, num_simulations) {
     #    purrr::map_dfr(1:num_simulations, function(sim_id) {
      #   individuals_genotypes <- initialize_individuals_pair(population, relationship_type, sim_id, loci_list)
       #  processed_genotypes <- log_function_time(process_individuals_genotypes, "process_individuals_genotypes", sim_id,individuals_genotypes, df_allelefreq, kinship_matrix)
        # return(processed_genotypes$result)
         #})
#     })
 #   fwrite(final_results, gsub("processed","v1_processed",output_file))
#  })
#}


int_results <- function(population, relationship_type, num_simulations){
                                     #, df_allelefreq, kinship_matrix, loci_list, output_file) {
   process_time <- system.time({
    f_results <- simulation_setup |>
      future_pmap_dfr(function(population, relationship_type, num_simulations){
        cs<-seq(from=1,to=num_simulations,by=1)

        foreach(i = cs, .combine="rbind") %dopar% {
          sim_id = cs[i]
          individuals_genotypes <- initialize_individuals_pair(population, relationship_type, sim_id, loci_list)

          processed_genotypes <- log_function_time(process_individuals_genotypes, "process_individuals_genotypes", sim_id,individuals_genotypes, df_allelefreq, kinship_matrix)
          return(processed_genotypes$result)
        }
      }, .options = furrr::furrr_options(seed = NULL))
    return(f_results)
  })

  log_message(paste("Processing completed in", process_time["elapsed"], "seconds."))

}

# Execute Simulation Setup and Processing (this is already too fast to benefit from being run in parallel)

setup_res <- log_function_time(generate_simulation_setup, "generate_simulation_setup", kinship_matrix, populations_list, n_sims_related, n_sims_unrelated)
simulation_setup <- setup_res$result

##########################################################################################################################################

if(use_remote_cluster==1){
  s_job<-slurm_apply(f=int_results, params=data.frame(simulation_setup), jobname = job_id,
                                          nodes = 1, cpus_per_node = 36, submit = TRUE,global_objects = c("initialize_individuals_pair","log_function_time","log_message"))
} else {
  job1<-slurm_apply(int_results, data.frame(simulation_setup), submit = FALSE,global_objects = c("initialize_individuals_pair","log_function_time","log_message"))
  s_job<-local_slurm_array(job1)
}

final_results<-get_slurm_out(s_job,"table",wait=TRUE)
fwrite(final_results, output_file,append=TRUE)
cleanup_files(s_job)

gc()
