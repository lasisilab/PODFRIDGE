# Load Required Libraries
library(dplyr)
library(furrr)
library(data.table)
library(future)
library(parallel)
library(doParallel)
library(stringr)

# Read Command-Line Arguments
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")
print(args)
population<-c("AfAm", "Cauc", "Hispanic", "Asian")

slurm_job_id <-  as.numeric(args[1])
job_name<-as.character(args[2])

# up the limit of memory available to future per core
options('future.globals.maxSize' = 1014*1024^2)

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

# Create output folder with SLURM job ID
output_dir <- file.path("output", paste0("simulation_", job_name))
dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, paste0("sim_processed_genotypes_",slurm_job_id,".csv"))
summary_output_file <- file.path(output_dir, "sim_summary_genotypes.csv")
timing_log_file <- file.path(output_dir, "timing_log.csv")

# Log the start of the process
log_message("Starting simulation setup and processing...")

###################################################################

# Load Allele Frequencies Data
log_message("Loading allele frequencies data...")
  allele_freq_time <- system.time({
  df_allelefreq <- fread(paste0(getwd(),"/data/df_allelefreq_combined.csv"))
  eval(parse(text=paste0("df_allelefreq <- df_allelefreq[df_allelefreq$population == \"",population[slurm_job_id],"\",]")))
  df_allelefreq[, allele := as.character(allele)]
})
log_message(paste("Loaded allele frequencies data in", allele_freq_time["elapsed"], "seconds."))
loci_list <- unique(df_allelefreq$marker)

# Load Core Loci Data
log_message("Loading core loci data...")
core_loci_time <- system.time({
  core_loci <- fread(paste0(getwd(),"/data/core_CODIS_loci.csv"))
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
rm(core_loci)

# Load Individuals Genotypes Data
input_dir <- file.path("data", "sims", paste0("simulation_script_", slurm_job_id,".out"))
input_dir <- file.path("data", paste0("simulation_", job_name))
individuals_genotypes <- fread(paste0(getwd(),"/",input_dir,"/processed_genotypes_",slurm_job_id,".csv"))
individuals_genotypes[,5:8] <- lapply(individuals_genotypes[,5:8],as.character)

# Define Kinship Matrix
kinship_matrix <- data.table(
  relationship_type = factor(c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins")),
  k0 = c(0, 1/4, 1/2, 7/8, 15/16),
  k1 = c(1, 1/2, 1/2, 1/8, 1/16),
  k2 = c(0, 1/4, 0, 0, 0)
)

# Functions
calculate_likelihood_ratio <- function(allele_frequency_data) {
  allele_frequency_data$LR<-ifelse(allele_frequency_data$shared_alleles == 0, allele_frequency_data$k0,NA)

  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1,allele_frequency_data$pA,NA)
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AA"|allele_frequency_data$genotype_match=="AA-AB",allele_frequency_data$Rxp*2,allele_frequency_data$Rxp)
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AC"|allele_frequency_data$genotype_match=="AB-AB",allele_frequency_data$Rxp*4,allele_frequency_data$Rxp)
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AB",(allele_frequency_data$Rxp*allele_frequency_data$pB)/(allele_frequency_data$pA+allele_frequency_data$pB),allele_frequency_data$Rxp)

  allele_frequency_data$LR<-ifelse(allele_frequency_data$shared_alleles == 0, allele_frequency_data$k0,NA)
  allele_frequency_data$LR<-ifelse(allele_frequency_data$shared_alleles == 1, allele_frequency_data$k0 + (allele_frequency_data$k1 / allele_frequency_data$Rxp),allele_frequency_data$LR)

                                   k0 = 0 #Why are we using these and not the values from the kinship matrix? Are they the same?
                                   k1 = 1
                                   k2 = 0

  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==2 & allele_frequency_data$genotype_match == "AA-AA",allele_frequency_data$pA,allele_frequency_data$Rxp)
  allele_frequency_data$Rxu<-ifelse(allele_frequency_data$shared_alleles==2  & allele_frequency_data$genotype_match == "AA-AA",allele_frequency_data$pA*allele_frequency_data$pA,NA)

  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==2 & allele_frequency_data$genotype_match == "AB-AB",(4*allele_frequency_data$pA*allele_frequency_data$pB)/(allele_frequency_data$pA+allele_frequency_data$pB),allele_frequency_data$Rxp)
  allele_frequency_data$Rxu<-ifelse(allele_frequency_data$shared_alleles==2  & allele_frequency_data$genotype_match == "AB-AB",2*allele_frequency_data$pA*allele_frequency_data$pB,allele_frequency_data$Rxu)

  allele_frequency_data$LR<-ifelse(allele_frequency_data$shared_alleles == 2, allele_frequency_data$k0 + (allele_frequency_data$k1 / allele_frequency_data$Rxp) + (allele_frequency_data$k2 / allele_frequency_data$Rxu),allele_frequency_data$LR)
  allele_frequency_data<-allele_frequency_data[, c("Rxp","Rxu","k0","k1","k2"):=NULL]

  return(allele_frequency_data)
}



kinship_calculation <- function(allele_frequency_data, kinship_matrix,df_allelefreq) {
  print("starting kinship_calculation")

  allele_frequency_data$alleles_ind1<- apply( allele_frequency_data[ , c('ind1_allele1','ind1_allele2') ] , 1 , paste , collapse = "_" )
  allele_frequency_data$alleles_ind2<- apply( allele_frequency_data[ , c('ind2_allele1','ind2_allele2') ] , 1 , paste , collapse = "_" )

  # Get shared alleles and their counts

  allele_frequency_data$shared_alleles_flag <- mapply(grepl, allele_frequency_data$ind1_allele1,allele_frequency_data$alleles_ind2)
  allele_frequency_data$shared_alleles1<-ifelse(allele_frequency_data$shared_alleles_flag,allele_frequency_data$ind1_allele1,NA)
  allele_frequency_data$shared_alleles_flag <- mapply(grepl, allele_frequency_data$ind1_allele2,allele_frequency_data$alleles_ind2)
  allele_frequency_data$shared_alleles1<-ifelse(allele_frequency_data$shared_alleles_flag & is.na(allele_frequency_data$shared_alleles1),allele_frequency_data$ind1_allele2,allele_frequency_data$shared_alleles1)

  allele_frequency_data$shared_alleles2<-ifelse(allele_frequency_data$shared_alleles_flag & !is.na(allele_frequency_data$shared_alleles1) & !allele_frequency_data$shared_alleles1==allele_frequency_data$ind1_allele2,allele_frequency_data$ind1_allele2,NA)

  allele_frequency_data$shared_alleles<-ifelse(!is.na(allele_frequency_data$shared_alleles1),1,0)
  allele_frequency_data$shared_alleles<-ifelse(!is.na(allele_frequency_data$shared_alleles2),2,allele_frequency_data$shared_alleles)

  allele_frequency_data<-allele_frequency_data[, c("shared_alleles_flag"):=NULL]

  #################################################################################################################################

 #Allele mapping
  
  allele_frequency_data$A<-ifelse(!is.na(allele_frequency_data$shared_alleles1),allele_frequency_data$shared_alleles1,allele_frequency_data$ind1_allele1)

  allele_frequency_data$B<-ifelse(!is.na(allele_frequency_data$shared_alleles2),allele_frequency_data$shared_alleles2,NA)
  allele_frequency_data$B<-ifelse(allele_frequency_data$shared_alleles==1 & !allele_frequency_data$ind1_allele1==allele_frequency_data$shared_alleles1,allele_frequency_data$ind1_allele1, allele_frequency_data$B)
  allele_frequency_data$B<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$ind1_allele1==allele_frequency_data$shared_alleles1 & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2,allele_frequency_data$ind1_allele2, allele_frequency_data$B)
  allele_frequency_data$B<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$ind1_allele1==allele_frequency_data$shared_alleles1 & allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & allele_frequency_data$ind2_allele1 == allele_frequency_data$shared_alleles1,allele_frequency_data$ind2_allele2, allele_frequency_data$B)
  allele_frequency_data$B<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$ind1_allele1==allele_frequency_data$shared_alleles1 & allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & !allele_frequency_data$ind2_allele1 == allele_frequency_data$shared_alleles1,allele_frequency_data$ind2_allele1, allele_frequency_data$B)

  allele_frequency_data$B<-ifelse(is.na(allele_frequency_data$shared_alleles1) & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2,allele_frequency_data$ind1_allele2, allele_frequency_data$B)
  allele_frequency_data$B<-ifelse(is.na(allele_frequency_data$shared_alleles1) & allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & !allele_frequency_data$ind2_allele1==allele_frequency_data$ind1_allele2,allele_frequency_data$ind2_allele1, allele_frequency_data$B)
  allele_frequency_data$B<-ifelse(allele_frequency_data$B==allele_frequency_data$A,allele_frequency_data$ind2_allele2, allele_frequency_data$B)
  allele_frequency_data$B<-ifelse(allele_frequency_data$B==allele_frequency_data$A,NA, allele_frequency_data$B)

  allele_frequency_data$C<-NA
  allele_frequency_data$C<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & !allele_frequency_data$ind2_allele1==allele_frequency_data$B,allele_frequency_data$ind2_allele1,allele_frequency_data$ind2_allele2)
  allele_frequency_data$C<-ifelse(allele_frequency_data$shared_alleles==1 & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & allele_frequency_data$ind2_allele1==allele_frequency_data$ind1_allele1,allele_frequency_data$ind2_allele2,allele_frequency_data$C)
  allele_frequency_data$C<-ifelse(allele_frequency_data$shared_alleles==1 & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & allele_frequency_data$ind2_allele1==allele_frequency_data$ind1_allele2,allele_frequency_data$ind2_allele2,allele_frequency_data$C)
  allele_frequency_data$C<-ifelse(allele_frequency_data$shared_alleles==1 & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & allele_frequency_data$ind2_allele2==allele_frequency_data$ind1_allele1,allele_frequency_data$ind2_allele1,allele_frequency_data$C)
  allele_frequency_data$C<-ifelse(allele_frequency_data$shared_alleles==1 & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & allele_frequency_data$ind2_allele2==allele_frequency_data$ind1_allele2,allele_frequency_data$ind2_allele1,allele_frequency_data$C)
  allele_frequency_data$C<-ifelse(allele_frequency_data$shared_alleles==0 & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2,allele_frequency_data$ind2_allele1,allele_frequency_data$C)
  allele_frequency_data$C<-ifelse(allele_frequency_data$shared_alleles==0 & allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2,allele_frequency_data$ind2_allele2,allele_frequency_data$C)
  allele_frequency_data$C<-ifelse(allele_frequency_data$C==allele_frequency_data$B,NA,allele_frequency_data$C)

  allele_frequency_data$D<-ifelse(allele_frequency_data$shared_alleles==0 & !allele_frequency_data$ind1_allele1==allele_frequency_data$ind1_allele2 & !allele_frequency_data$ind2_allele1==allele_frequency_data$ind2_allele2,allele_frequency_data$ind2_allele2,NA)

 # Generate genotype strings
 allele_frequency_data$labeled_alleles_ind1a<-ifelse(allele_frequency_data$ind1_allele1==allele_frequency_data$A,"A",NA)
 allele_frequency_data$labeled_alleles_ind1a<-ifelse(is.na(allele_frequency_data$labeled_alleles_ind1a) & allele_frequency_data$ind1_allele1==allele_frequency_data$B,"B",allele_frequency_data$labeled_alleles_ind1a)

 allele_frequency_data$labeled_alleles_ind1b<-ifelse(allele_frequency_data$ind1_allele2==allele_frequency_data$A,"A",NA)
 allele_frequency_data$labeled_alleles_ind1b<-ifelse(is.na(allele_frequency_data$labeled_alleles_ind1b) & allele_frequency_data$ind1_allele2==allele_frequency_data$B,"B",allele_frequency_data$labeled_alleles_ind1b)

 allele_frequency_data$labeled_alleles_ind2a<-ifelse(allele_frequency_data$ind2_allele1==allele_frequency_data$A,"A",NA)
 allele_frequency_data$labeled_alleles_ind2a<-ifelse(is.na(allele_frequency_data$labeled_alleles_ind2a) & allele_frequency_data$ind2_allele1==allele_frequency_data$B,"B",allele_frequency_data$labeled_alleles_ind2a)
 allele_frequency_data$labeled_alleles_ind2a<-ifelse(is.na(allele_frequency_data$labeled_alleles_ind2a) & allele_frequency_data$ind2_allele1==allele_frequency_data$C,"C",allele_frequency_data$labeled_alleles_ind2a)

 allele_frequency_data$labeled_alleles_ind2b<-ifelse(allele_frequency_data$ind2_allele2==allele_frequency_data$A,"A",NA)
 allele_frequency_data$labeled_alleles_ind2b<-ifelse(is.na(allele_frequency_data$labeled_alleles_ind2b) & allele_frequency_data$ind2_allele2==allele_frequency_data$B,"B",allele_frequency_data$labeled_alleles_ind2b)
 allele_frequency_data$labeled_alleles_ind2b<-ifelse(is.na(allele_frequency_data$labeled_alleles_ind2b) & allele_frequency_data$ind2_allele2==allele_frequency_data$C,"C",allele_frequency_data$labeled_alleles_ind2b)
 allele_frequency_data$labeled_alleles_ind2b<-ifelse(is.na(allele_frequency_data$labeled_alleles_ind2b) & allele_frequency_data$ind2_allele2==allele_frequency_data$D,"D",allele_frequency_data$labeled_alleles_ind2b)

 allele_frequency_data$which_1st <- with(allele_frequency_data,labeled_alleles_ind1a<labeled_alleles_ind1b)
 allele_frequency_data$genotype_ind1 <- with(allele_frequency_data,  ifelse(which_1st,paste(labeled_alleles_ind1a,labeled_alleles_ind1b,sep=""),paste(labeled_alleles_ind1b,labeled_alleles_ind1a,sep="")))

 allele_frequency_data$which_1st <- with(allele_frequency_data,labeled_alleles_ind2a<labeled_alleles_ind2b)
 allele_frequency_data$genotype_ind2 <- with(allele_frequency_data,
                                             ifelse(which_1st,paste(labeled_alleles_ind2a,labeled_alleles_ind2b,sep=""),paste(labeled_alleles_ind2b,labeled_alleles_ind2a,sep="")))

 allele_frequency_data$genotype_match <- paste(allele_frequency_data$genotype_ind1, allele_frequency_data$genotype_ind2, sep = "-")
 allele_frequency_data<-allele_frequency_data[, c("labeled_alleles_ind2a","labeled_alleles_ind2b","labeled_alleles_ind1a","labeled_alleles_ind1b"):=NULL]

 ########################################################################################################################################################################

  allele_frequency_data<-left_join(allele_frequency_data,kinship_matrix)
  allele_frequency_data$marker<-allele_frequency_data$locus

  t1<-df_allelefreq
  names(t1)<-c("A","locus","pA","population")
  allele_frequency_data<-left_join(allele_frequency_data,t1)

  t1<-df_allelefreq
  names(t1)<-c("B","locus","pB","population")
  allele_frequency_data<-left_join(allele_frequency_data,t1)
  rm(t1)

  allele_frequency_data$relationship_known<-allele_frequency_data$relationship_type
  allele_frequency_data$relationship_tested<-allele_frequency_data$relationship_type
  
print("kinship_calculations")
  kinship_calculations <- calculate_likelihood_ratio(allele_frequency_data)
  return(kinship_calculations)
}

calculate_combined_lrs <- function(final_results, loci_lists) {
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

log_message("Processing individuals genotypes...")

processing_time <- system.time({
  processed_genotypes <- kinship_calculation(individuals_genotypes, kinship_matrix,df_allelefreq)
})

log_message(paste("Processed individuals genotypes in", processing_time["elapsed"], "seconds."))
processed_genotypes <- as.data.table(processed_genotypes)

# Calculate combined likelihood ratios - all data from each population needs to be together at this point
log_message("Calculating combined likelihood ratios...")
  combined_lrs_time <- system.time({
  combined_lrs <- calculate_combined_lrs(processed_genotypes, loci_lists)
})
  rm(loci_lists)
log_message(paste("Calculated combined likelihood ratios in", combined_lrs_time["elapsed"], "seconds."))

# Save results to CSV
log_message("Saving results to CSV files...")
fwrite(processed_genotypes, output_file)
fwrite(combined_lrs, summary_output_file,append=TRUE)

# Save timing log to CSV
timing_log_df <- as.data.frame(rbind(c("combined_lrs",combined_lrs_time), c("genotype_processing",processing_time)))

fwrite(timing_log_df, timing_log_file)

log_message("Simulation processing completed.")

