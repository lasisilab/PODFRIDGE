# Load Required Libraries
suppressMessages(suppressWarnings({
  library(tidyverse)
  library(furrr)
  library(data.table)
  library(future)
  library(parallel)
  library(doParallel)
  library(benchmarkme)
}))

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
  plan(cluster, workers = cl)
} else { #if sending to remote cluster(s)
  plan(cluster, workers = ("clustername1")) #use syntax "server.remote.org" in workers if using an online cluster
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
args <- commandArgs(trailingOnly = TRUE)
#args<-c("100","200","100200_runthrough_focal")
n_sims_related <- as.numeric(args[1])
n_sims_unrelated <- as.numeric(args[2])
job_id <- as.character(args[3])

# Create output folder with job ID
output_dir <- file.path("data", "sims", paste0("simulation_script_", slurm_job_id,".out"))
dir.create(output_dir, recursive = TRUE)
output_file <- "processed_genotypes.csv"
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
  df_allelefreq<-df_allelefreq %>%
    tidyr::complete(population,allele,marker)

})
log_message(paste("Loaded allele frequencies data in", allele_freq_time["elapsed"], "seconds."))

test1<-df_allelefreq[df_allelefreq$population=="AfAm",]
test2<-df_allelefreq[df_allelefreq$population=="Asian",]

#Replace 0 frequencies of absent allele and population combinations - assign the rare allele frequency 5/(2*1036)
df_allelefreq <- replace(df_allelefreq, df_allelefreq==0, 5/(2*1036))
df_allelefreq <- replace(df_allelefreq, is.na(df_allelefreq), 5/(2*1036))

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
num_1s<-round(n_sims_related/length(populations_list),0)

print("Defining functions")

initialize_individuals <- function(populations_list, num_related, loci_list) {
   num_loci <- length(loci_list) #Total number of loci- around 1000
   individuals_genotypes2<-foreach(i = 1:length(populations_list),.export=c("loci_list","num_loci","num_1s"),.combine="rbind") %dopar% {
    individuals_genotypes <- data.table(
    sim_id = rep(1:num_related, each=num_loci),#repeat 1:length simid, each num_loci length.
    locus = rep(loci_list, num_related)) #repeat each once for n sim_id
    individuals_genotypes$population<-populations_list[i]
    return(individuals_genotypes)
  }
  return(individuals_genotypes2)
}

###############################################################################################################################

simulate_genotypes1 <- function(row, df_allelefreq, kinship_matrix,related) {
  population <- row$population
  locus <- row$locus
  allele_freqs <- df_allelefreq[which(df_allelefreq$population == population & df_allelefreq$marker == locus),]
  alleles <- allele_freqs$allele
  frequencies <- allele_freqs$frequency
  frequencies <- round(frequencies / sum(frequencies), 6)
  valid_indices <- frequencies > 0
  alleles <- alleles[valid_indices]
  frequencies <- frequencies[valid_indices]
  ind1_alleles <- sample(alleles, size = 2, replace = TRUE, prob = frequencies)
  subtab<-row
  if(related==1){
    kinship_matrix<-kinship_matrix[!kinship_matrix$relationship_type=="unrelated",]
    subtab<-cbind(subtab, rep(row.names(subtab), each = length(kinship_matrix$relationship_type)))
  } else {
    kinship_matrix<-kinship_matrix[kinship_matrix$relationship_type=="unrelated",]
    subtab$relationship_type<-"unrelated"
  }
  subtab$relationship_type<-kinship_matrix$relationship_type
  subtab$ind1_allele1 <- ind1_alleles[1]
  subtab$ind1_allele2 <- ind1_alleles[2]
 return(subtab)
}

simulate_genotypes2 <- function(row, df_allelefreq, kinship_matrix) {

  population <- row$population
  locus <- row$locus
  relationship <- as.factor(row$relationship_type)

  allele_freqs <- df_allelefreq[which(df_allelefreq$population == population & df_allelefreq$marker == locus),]
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

  row$ind2_allele1 <- ind2_alleles[1]
  row$ind2_allele2 <- ind2_alleles[2]
  return(row)
}

process_simulation_setup <- function(df_allelefreq, kinship_matrix, loci_list, populations_list,num_related,num_unrelated) {
  process_time <- system.time({
    #Start by generating the 'individual 1's for each sim per population
     individuals_genotypes <- initialize_individuals(populations_list, num_related, loci_list)
     individuals_genotypes$relationship_type<-"related" #These are placeholder fields which are later overwritten
     individuals_genotypes2 <- initialize_individuals(populations_list, num_unrelated, loci_list)
     individuals_genotypes2$relationship_type<-"unrelated"
     n1<-names(individuals_genotypes)
     ##############################################################################################################

     #Add 'individual 1' genotypes row by row
     final_resultsa <- individuals_genotypes |>
       furrr::future_pmap(function(...) {
         processed_genotypes<- log_function_time(simulate_genotypes1, "simulate_genotypes1", list(...), df_allelefreq, kinship_matrix,1)
         return(processed_genotypes$result)
             }, seed = TRUE) |>
       data.table::rbindlist()

      final_resultsb <- individuals_genotypes2 |>
        furrr::future_pmap(function(...) {
          processed_genotypes<- log_function_time(simulate_genotypes1, "simulate_genotypes1", list(...), df_allelefreq, kinship_matrix,0)
          return(processed_genotypes$result)
        }, seed = TRUE) |>
        data.table::rbindlist()

      final_resultsa<-final_resultsa[,c(1:3,6:8)]
      final_resultsb<-final_resultsb[,c(1:4,6:7)]

      names(final_resultsa)<-names(final_resultsb)
    final_results<-rbind(final_resultsa,final_resultsb)
    return(final_results)
  })
  log_message(paste("Processing completed in", process_time["elapsed"], "seconds."))

}

process_simulation_setup2 <- function(individuals_genotypes,df_allelefreq, kinship_matrix,output_file) {
  process_time <- system.time({

  final_resultsc <- individuals_genotypes |>
    furrr::future_pmap(function(...) {
      processed_genotypes<- log_function_time(simulate_genotypes2, "simulate_genotypes2", list(...), df_allelefreq, kinship_matrix)
      return(processed_genotypes$result)
    }, seed = TRUE) |>
    data.table::rbindlist()
  final_resultsc <-final_resultsc[, c("seed"):=NULL]
   #  return(final_resultsc)
  #  fwrite(final_resultsc, output_file,append=TRUE)
  })
  log_message(paste("Processing completed in", process_time["elapsed"], "seconds."))
  return(final_resultsc)
}


#########################################################################################################################

print("Generating simulations")
proc_res1 <- log_function_time(process_simulation_setup, "process_simulation_setup",df_allelefreq, kinship_matrix,loci_list,populations_list,n_sims_related,n_sims_unrelated)
timing_log <- proc_res1$timing_log

proc_res1<-data.table(proc_res1$result)

proc_res2 <- log_function_time(process_simulation_setup2, "process_simulation_setup2", proc_res1, df_allelefreq, kinship_matrix, output_file)
timing_log <- proc_res2$timing_log
#proc_res2<-proc_res2$result
final_resultsc <-as.data.frame(proc_res2$result)
write.table(final_resultsc,gzfile(paste0(output_dir,"/",output_file,".gz")),append=FALSE)
#save(final_resultsc,file=output_file, compress=T)
#write.csv(final_resultsc,"final_resultsc.csv")


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
#gc()
