#Based on lr.R
#This script expands population and relationship options to create additional LR columns rather than rows
#LR assumes known and tested populations are the same, and that known and tested relationships are accurate
#LR columns specifying a different predicted population assume that the predicted relationship is accurate
#LR columns specifying a different predicted relationship assume that the predicted population is accurate
#Combined LRs are based on the original LR values
#Combined LR fields with a '_tested' suffix show alternative LRs based on the tested (not known) relationships

# Load Required Libraries
library(dplyr)
library(furrr)
library(data.table)
library(future)
library(parallel)
library(doParallel)
library(stringr)

# Read Command-Line Arguments
#args<-c("3","simulation_100200_runthrough_focal")
args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")

slurm_job_id <-  as.numeric(args[1])
job_name<-as.character(args[2])
print(slurm_job_id)
print(job_name)
t<-getwd()
# up the limit of memory available to future per core
options('future.globals.maxSize' = 1014*1024^2)

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}

print(getwd())

# Create output folder with SLURM job ID
t<-getwd()

input_dir <- file.path("data", "sims", paste0("simulation_script_", slurm_job_id,".out"))
dir.create(input_dir, recursive = TRUE)
output_dir <- input_dir
#output_file <- file.path(output_dir, paste0("sim_processed_genotypes.csv"))
output_file <- "sim_processed_genotypes.csv"
output_file2 <- "sim_combined_genotypes.csv"
timing_log_file <- file.path(output_dir, "timing_log.csv")

# Log the start of the process
log_message("Starting LR processing for known relationships")


###################################################################

# Load Allele Frequencies Data
log_message("Loading allele frequencies data...")
print(getwd())
allele_freq_time <- system.time({
  df_allelefreq <- fread(paste0(getwd(),"/data/df_allelefreq_combined.csv"))
  df_allelefreq$frequency <- ifelse(df_allelefreq$frequency==0,5/(2*1036),df_allelefreq$frequency)
  df_allelefreq <-df_allelefreq[!df_allelefreq$population=="all",]

  #  eval(parse(text=paste0("df_allelefreq <- df_allelefreq[df_allelefreq$population == \"",population[slurm_job_id],"\",]")))
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
#individuals_genotypes <- fread(paste0(input_dir,"/processed_genotypes.csv"))
#individuals_genotypes <- load(paste0(input_dir,"/processed_genotypes.Rbin"))
individuals_genotypes <- read.table(gzfile(paste0(input_dir,"/processed_genotypes.csv.gz")),row.names=1)
individuals_genotypes <-individuals_genotypes %>%
  select(-one_of("seed")) %>%
  as.data.table()

individuals_genotypes[,5:8] <- lapply(individuals_genotypes[,5:8],as.character)
names(individuals_genotypes)[3]<-"population_known"

# Define Kinship Matrix
kinship_matrix <- data.table(
  relationship_type = factor(c("parent_child", "full_siblings", "half_siblings", "cousins", "second_cousins","unrelated")),
  k0 = c(0, 1/4, 1/2, 7/8, 15/16, 1),
  k1 = c(1, 1/2, 1/2, 1/8, 1/16, 0),
  k2 = c(0, 1/4, 0, 0, 0, 0)
)

#################################################################################

# Functions
calculate_likelihood_ratio <- function(allele_frequency_data,population) {

  for(p in 1:length(population)){
    eval(parse(text=paste0("allele_frequency_data$LR_",population[p],"<-ifelse(allele_frequency_data$shared_alleles == 0, allele_frequency_data$k0,NA)")))
    eval(parse(text=paste0("allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1,allele_frequency_data$pA_",population[p],",NA)")))
    eval(parse(text=paste0("allele_frequency_data$pA <- allele_frequency_data$pA_",population[p])))
    eval(parse(text=paste0("allele_frequency_data$pB <- allele_frequency_data$pB_",population[p])))
    allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1,allele_frequency_data$pA,NA)
    allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AA"|allele_frequency_data$genotype_match=="AA-AB",allele_frequency_data$Rxp*2,allele_frequency_data$Rxp)
    allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AC"|allele_frequency_data$genotype_match=="AB-AB",allele_frequency_data$Rxp*4,allele_frequency_data$Rxp)
    allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AB",(allele_frequency_data$Rxp*allele_frequency_data$pB)/(allele_frequency_data$pA+allele_frequency_data$pB),allele_frequency_data$Rxp)
    eval(parse(text=paste0("allele_frequency_data$LR_",population[p],"<-ifelse(allele_frequency_data$shared_alleles == 1, allele_frequency_data$k0 + (allele_frequency_data$k1 / allele_frequency_data$Rxp),allele_frequency_data$LR_",population[p],")")))

    allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==2 & allele_frequency_data$genotype_match == "AA-AA",allele_frequency_data$pA,allele_frequency_data$Rxp)
    allele_frequency_data$Rxu<-ifelse(allele_frequency_data$shared_alleles==2  & allele_frequency_data$genotype_match == "AA-AA",allele_frequency_data$pA*allele_frequency_data$pA,NA)

    allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==2 & allele_frequency_data$genotype_match == "AB-AB",(4*allele_frequency_data$pA*allele_frequency_data$pB)/(allele_frequency_data$pA+allele_frequency_data$pB),allele_frequency_data$Rxp)
    allele_frequency_data$Rxu<-ifelse(allele_frequency_data$shared_alleles==2  & allele_frequency_data$genotype_match == "AB-AB",2*allele_frequency_data$pA*allele_frequency_data$pB,allele_frequency_data$Rxu)

    eval(parse(text=paste0("allele_frequency_data$LR_",population[p],"<-ifelse(allele_frequency_data$shared_alleles == 2,
      allele_frequency_data$k0 + (allele_frequency_data$k1 / allele_frequency_data$Rxp) + (allele_frequency_data$k2 / allele_frequency_data$Rxu),
      allele_frequency_data$LR_",population[p],")")))

  }


  allele_frequency_data<-allele_frequency_data[, c("Rxp","Rxu","k0","k1","k2"):=NULL]
  allele_frequency_data<-allele_frequency_data[, c("pA","pB","shared_alleles1","shared_alleles2","which_1stb"):=NULL]
  allele_frequency_data<-allele_frequency_data[, c("pA_AfAm","pB_AfAm","pA_Asian","pB_Asian","pA_Cauc","pB_Cauc","pA_Hispanic","pB_Hispanic"):=NULL]
  return(allele_frequency_data)
}

calculate_likelihood_ratio2 <- function(allele_frequency_data) {
  allele_frequency_data$LR_relationship_tested<-ifelse(allele_frequency_data$shared_alleles == 0, allele_frequency_data$k0,NA)
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1,allele_frequency_data$pA,NA)
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AA"|allele_frequency_data$genotype_match=="AA-AB",allele_frequency_data$Rxp*2,allele_frequency_data$Rxp)
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AC"|allele_frequency_data$genotype_match=="AB-AB",allele_frequency_data$Rxp*4,allele_frequency_data$Rxp)
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==1 & allele_frequency_data$genotype_match=="AB-AB",(allele_frequency_data$Rxp*allele_frequency_data$pB)/(allele_frequency_data$pA+allele_frequency_data$pB),allele_frequency_data$Rxp)
  allele_frequency_data$LR_relationship_tested<-ifelse(allele_frequency_data$shared_alleles == 1, allele_frequency_data$k0 + (allele_frequency_data$k1 / allele_frequency_data$Rxp),allele_frequency_data$LR_relationship_tested)

  k0 = 0 #Why are we using these and not the values from the kinship matrix? Are they the same?
  k1 = 1
  k2 = 0
  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==2 & allele_frequency_data$genotype_match == "AA-AA",allele_frequency_data$pA,allele_frequency_data$Rxp)
  allele_frequency_data$Rxu<-ifelse(allele_frequency_data$shared_alleles==2  & allele_frequency_data$genotype_match == "AA-AA",allele_frequency_data$pA*allele_frequency_data$pA,NA)

  allele_frequency_data$Rxp<-ifelse(allele_frequency_data$shared_alleles==2 & allele_frequency_data$genotype_match == "AB-AB",(4*allele_frequency_data$pA*allele_frequency_data$pB)/(allele_frequency_data$pA+allele_frequency_data$pB),allele_frequency_data$Rxp)
  allele_frequency_data$Rxu<-ifelse(allele_frequency_data$shared_alleles==2  & allele_frequency_data$genotype_match == "AB-AB",2*allele_frequency_data$pA*allele_frequency_data$pB,allele_frequency_data$Rxu)

  allele_frequency_data$LR_relationship_tested<-ifelse(allele_frequency_data$shared_alleles == 2,
  allele_frequency_data$k0 + (allele_frequency_data$k1 / allele_frequency_data$Rxp) + (allele_frequency_data$k2 / allele_frequency_data$Rxu),
  allele_frequency_data$LR_relationship_tested)

  allele_frequency_data<-allele_frequency_data[, c("Rxp","Rxu","k0","k1","k2"):=NULL]
  allele_frequency_data<-allele_frequency_data[, c("pA","pB"):=NULL]
  return(allele_frequency_data)
}


kinship_calculation <- function(allele_frequency_data, kinship_matrix,df_allelefreq) {
  print("starting kinship_calculation")

  allele_frequency_data$alleles_ind1<- apply( allele_frequency_data[ , c('ind1_allele1','ind1_allele2') ] , 1 , paste , collapse = "_" )
  allele_frequency_data$alleles_ind2<- apply( allele_frequency_data[ , c('ind2_allele1','ind2_allele2') ] , 1 , paste , collapse = "_" )

  # Get shared alleles and their counts

  allele_frequency_data$shared_alleles1<-NA
  allele_frequency_data$shared_alleles_flag <- ifelse(allele_frequency_data$ind1_allele1==allele_frequency_data$ind2_allele1|allele_frequency_data$ind1_allele1==allele_frequency_data$ind2_allele2,1,0)
  allele_frequency_data$shared_alleles1<-ifelse(allele_frequency_data$shared_alleles_flag==1,allele_frequency_data$ind1_allele1,NA)
  allele_frequency_data$shared_alleles_flag <- ifelse(allele_frequency_data$ind1_allele2==allele_frequency_data$ind2_allele1|allele_frequency_data$ind1_allele2==allele_frequency_data$ind2_allele2,1,0)
  allele_frequency_data$shared_alleles1<-ifelse(allele_frequency_data$shared_alleles_flag==1 & is.na(allele_frequency_data$shared_alleles1),
                                                allele_frequency_data$ind1_allele2,allele_frequency_data$shared_alleles1)
  allele_frequency_data$shared_alleles2<-ifelse(allele_frequency_data$shared_alleles_flag==1 & !is.na(allele_frequency_data$shared_alleles1) &
                                                  !allele_frequency_data$shared_alleles1==allele_frequency_data$ind1_allele2,
                                                allele_frequency_data$ind1_allele2,NA)

  allele_frequency_data$shared_alleles<-ifelse(!is.na(allele_frequency_data$shared_alleles1),1,0)
  allele_frequency_data$shared_alleles<-ifelse(!is.na(allele_frequency_data$shared_alleles2),2,allele_frequency_data$shared_alleles)

  allele_frequency_data<-allele_frequency_data[, c("shared_alleles_flag"):=NULL]

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

  allele_frequency_data$which_1stb <- with(allele_frequency_data,labeled_alleles_ind2a<labeled_alleles_ind2b)
  allele_frequency_data$genotype_ind2 <- with(allele_frequency_data,
                                              ifelse(which_1stb,paste(labeled_alleles_ind2a,labeled_alleles_ind2b,sep=""),paste(labeled_alleles_ind2b,labeled_alleles_ind2a,sep="")))

  allele_frequency_data$genotype_match <- paste(allele_frequency_data$genotype_ind1, allele_frequency_data$genotype_ind2, sep = "-")
  allele_frequency_data<-allele_frequency_data[, c("labeled_alleles_ind2a","labeled_alleles_ind2b","labeled_alleles_ind1a","labeled_alleles_ind1b","seed"):=NULL]


  #- for each row we need the marker that matches the locus, extract the frequencies from the whole table for the A and B alleles

  #################################################################
  #PART A-
  #Assuming our relationships tested and known are the same-
  #Split results by population- create each combination of predicted and known population-
  #This cannot be done in long format as this is too large so create a wide format reference table to join on

  allele_frequency_data<-left_join(allele_frequency_data,kinship_matrix)
  allele_frequency_data$marker<-allele_frequency_data$locus
  population<-c("AfAm", "Cauc", "Hispanic", "Asian")
  if(exists("df3")){
    rm(df3)
  }

  df2<-df_allelefreq[df_allelefreq$population=="AfAm",]
  names(df2)[3]<-"frequency_AfAm"
  df2<-df2 %>%
    dplyr::select(-'population')
  df3<-df_allelefreq[df_allelefreq$population=="Asian",]
  names(df3)[3]<-"frequency_Asian"
  df3<-df3 %>%
    dplyr::select(-'population')
  df4<-df_allelefreq[df_allelefreq$population=="Cauc",]
  names(df4)[3]<-"frequency_Cauc"
  df4<-df4 %>%
    dplyr::select(-'population')
  df5<-df_allelefreq[df_allelefreq$population=="Hispanic",]
  names(df5)[3]<-"frequency_Hispanic"
  df5<-df5 %>%
    dplyr::select(-'population')
  t1<-left_join(df2,df3)
  t1<-left_join(t1,df4)
  t1<-left_join(t1,df5)

  names(t1)<-gsub("frequency","pA",names(t1))
  names(t1)[1:2]<-c("A","locus")
  allele_frequency_data<-left_join(allele_frequency_data,t1)

  t2<-t1
  names(t2)<-gsub("frequency","pB",names(t2))
  names(t2)[1:2]<-c("B","locus")
  allele_frequency_data<-left_join(allele_frequency_data,t2)
#  summary(allele_frequency_data)
  rm(t1)
  rm(t2)
  rm(df2)
  rm(df3)

  allele_frequency_data <- calculate_likelihood_ratio(allele_frequency_data,population)
  allele_frequency_data$LR <- NA
  allele_frequency_data$LR <- ifelse(allele_frequency_data$population_known=="AfAm",allele_frequency_data$LR_AfAm,allele_frequency_data$LR)
  allele_frequency_data$LR <- ifelse(allele_frequency_data$population_known=="Asian",allele_frequency_data$LR_Asian,allele_frequency_data$LR)
  allele_frequency_data$LR <- ifelse(allele_frequency_data$population_known=="Cauc",allele_frequency_data$LR_Cauc,allele_frequency_data$LR)
  allele_frequency_data$LR <- ifelse(allele_frequency_data$population_known=="Hispanic",allele_frequency_data$LR_Hispanic,allele_frequency_data$LR)

  allele_frequency_data$LR_incorrect_population<-NA
  allele_frequency_data$LR_incorrect_population <- ifelse(allele_frequency_data$population_known=="AfAm",(allele_frequency_data$LR_Hispanic*allele_frequency_data$LR_Asian*allele_frequency_data$LR_Cauc),allele_frequency_data$LR_incorrect_population)
  allele_frequency_data$LR_incorrect_population <- ifelse(allele_frequency_data$population_known=="Asian",prod(allele_frequency_data$LR_AfAm*allele_frequency_data$LR_Hispanic*allele_frequency_data$LR_Cauc),allele_frequency_data$LR_incorrect_population)
  allele_frequency_data$LR_incorrect_population <- ifelse(allele_frequency_data$population_known=="Cauc",prod(allele_frequency_data$LR_AfAm*allele_frequency_data$LR_Asian*allele_frequency_data$LR_Hispanic),allele_frequency_data$LR_incorrect_population)
  allele_frequency_data$LR_incorrect_population <- ifelse(allele_frequency_data$population_known=="Hispanic",prod(allele_frequency_data$LR_AfAm*allele_frequency_data$LR_Asian*allele_frequency_data$LR_Cauc),allele_frequency_data$LR_incorrect_population)


  #################################################################
  #PART B-
  #Assuming our populations tested and known are the same
  #Split results by tested relationship- create each combination of predicted and known relationship
  #This creates a new row for each combination of relationship tested vs known

  names(allele_frequency_data)[4]<-"relationship_known"

  relationship_type<-kinship_matrix$relationship_type
  length_r<-nrow(kinship_matrix)
  length_t<-nrow(allele_frequency_data)

  allele_frequency_data<-sapply(allele_frequency_data, rep.int, times=length_r)
  allele_frequency_data<-as.data.table(allele_frequency_data)
  relationship_tested<-rep(kinship_matrix$relationship_type,each=length_t)
  allele_frequency_data<-cbind(allele_frequency_data,relationship_tested)

  names(kinship_matrix)[1]<-"relationship_tested"
  allele_frequency_data<-left_join(allele_frequency_data,kinship_matrix,by="relationship_tested")

  #Add on the pA and pB assuming the 'correct' population



  t1<-df_allelefreq
  names(t1)<-c("A","marker","pA","population_known")
  allele_frequency_data<-left_join(allele_frequency_data,t1)

  t2<-df_allelefreq
  names(t2)<-c("B","marker","pB","population_known")
  allele_frequency_data<-left_join(allele_frequency_data,t2)

  rm(t1)
  rm(t2)

  kinship_calculations <- calculate_likelihood_ratio2(allele_frequency_data)

  return(kinship_calculations)
}


calculate_combined_lrs <- function(final_results, loci_lists) {
  final_results[sapply(final_results, is.infinite)] <- 0

  #The 'LR' column assumes the correct population has been tested
  print("Run calculate_combined_lrs")

 final_results$locus<-as.character(final_results$locus)
 final_results$LR<-as.numeric(final_results$LR)
 final_results$LR_relationship_tested<-as.numeric(final_results$LR_relationship_tested)
 final_results$LR_incorrect_population<-as.numeric(final_results$LR_incorrect_population)

  combined_lrs <- final_results[, .(
    core_13 = prod(LR[locus %in% loci_lists$core_13], na.rm = TRUE),
    identifiler_15 = prod(LR[locus %in% loci_lists$identifiler_15], na.rm = TRUE),
    expanded_20 = prod(LR[locus %in% loci_lists$expanded_20], na.rm = TRUE),
    supplementary = prod(LR[locus %in% loci_lists$supplementary], na.rm = TRUE),
    autosomal_29 = prod(LR[locus %in% loci_lists$autosomal_29], na.rm = TRUE)
  ), by = .(population_known, relationship_known, relationship_tested, sim_id)]
  print(head(combined_lrs))
  print(2)
  combined_lrs <- melt(combined_lrs,
                       id.vars = c("population_known", "relationship_known", "relationship_tested", "sim_id"),
                       measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                       variable.name = "loci_set", value.name = "LR")

  combined_lrs2 <- final_results[, .(
  core_13 = prod(LR_relationship_tested[locus %in% loci_lists$core_13], na.rm = TRUE),
  identifiler_15 = prod(LR_relationship_tested[locus %in% loci_lists$identifiler_15], na.rm = TRUE),
  expanded_20 = prod(LR_relationship_tested[locus %in% loci_lists$expanded_20], na.rm = TRUE),
  supplementary = prod(LR_relationship_tested[locus %in% loci_lists$supplementary], na.rm = TRUE),
  autosomal_29 = prod(LR_relationship_tested[locus %in% loci_lists$autosomal_29], na.rm = TRUE)
  ), by = .(population_known, relationship_known, relationship_tested, sim_id)]


  final_results$LR_AfAm<-as.numeric(final_results$LR_AfAm)
  final_results$LR_Asian<-as.numeric(final_results$LR_Asian)
  final_results$LR_Cauc<-as.numeric(final_results$LR_Cauc)
  final_results$LR_Hispanic<-as.numeric(final_results$LR_Hispanic)

  combined_lrs3 <- final_results[, .(
    core_13 = prod(LR_AfAm[locus %in% loci_lists$core_13], na.rm = TRUE),
    identifiler_15 = prod(LR_AfAm[locus %in% loci_lists$identifiler_15], na.rm = TRUE),
    expanded_20 = prod(LR_AfAm[locus %in% loci_lists$expanded_20], na.rm = TRUE),
    supplementary = prod(LR_AfAm[locus %in% loci_lists$supplementary], na.rm = TRUE),
    autosomal_29 = prod(LR_AfAm[locus %in% loci_lists$autosomal_29], na.rm = TRUE)
  ), by = .(population_known, relationship_known, relationship_tested, sim_id)]

  combined_lrs4 <- final_results[, .(
    core_13 = prod(LR_Asian[locus %in% loci_lists$core_13], na.rm = TRUE),
    identifiler_15 = prod(LR_Asian[locus %in% loci_lists$identifiler_15], na.rm = TRUE),
    expanded_20 = prod(LR_Asian[locus %in% loci_lists$expanded_20], na.rm = TRUE),
    supplementary = prod(LR_Asian[locus %in% loci_lists$supplementary], na.rm = TRUE),
    autosomal_29 = prod(LR_Asian[locus %in% loci_lists$autosomal_29], na.rm = TRUE)
  ), by = .(population_known, relationship_known, relationship_tested, sim_id)]

  combined_lrs5 <- final_results[, .(
    core_13 = prod(LR_Cauc[locus %in% loci_lists$core_13], na.rm = TRUE),
    identifiler_15 = prod(LR_Cauc[locus %in% loci_lists$identifiler_15], na.rm = TRUE),
    expanded_20 = prod(LR_Cauc[locus %in% loci_lists$expanded_20], na.rm = TRUE),
    supplementary = prod(LR_Cauc[locus %in% loci_lists$supplementary], na.rm = TRUE),
    autosomal_29 = prod(LR_Cauc[locus %in% loci_lists$autosomal_29], na.rm = TRUE)
  ), by = .(population_known, relationship_known, relationship_tested, sim_id)]

  combined_lrs6 <- final_results[, .(
    core_13 = prod(LR_Hispanic[locus %in% loci_lists$core_13], na.rm = TRUE),
    identifiler_15 = prod(LR_Hispanic[locus %in% loci_lists$identifiler_15], na.rm = TRUE),
    expanded_20 = prod(LR_Hispanic[locus %in% loci_lists$expanded_20], na.rm = TRUE),
    supplementary = prod(LR_Hispanic[locus %in% loci_lists$supplementary], na.rm = TRUE),
    autosomal_29 = prod(LR_Hispanic[locus %in% loci_lists$autosomal_29], na.rm = TRUE)
  ), by = .(population_known, relationship_known, relationship_tested, sim_id)]

  combined_lrs2 <- melt(combined_lrs2,
                       id.vars = c("population_known", "relationship_known", "relationship_tested", "sim_id"),
                       measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                       variable.name = "loci_set", value.name = "LR_relationship_tested")

  combined_lrs3 <- melt(combined_lrs3,
                        id.vars = c("population_known", "relationship_known", "relationship_tested", "sim_id"),
                        measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                        variable.name = "loci_set", value.name = "LR_Afam")

  combined_lrs4 <- melt(combined_lrs4,
                        id.vars = c("population_known", "relationship_known", "relationship_tested", "sim_id"),
                        measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                        variable.name = "loci_set", value.name = "LR_Asian")

  combined_lrs5 <- melt(combined_lrs5,
                        id.vars = c("population_known", "relationship_known", "relationship_tested", "sim_id"),
                        measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                        variable.name = "loci_set", value.name = "LR_Cauc")

  combined_lrs6 <- melt(combined_lrs6,
                        id.vars = c("population_known", "relationship_known", "relationship_tested", "sim_id"),
                        measure.vars = c("core_13", "identifiler_15", "expanded_20", "supplementary", "autosomal_29"),
                        variable.name = "loci_set", value.name = "LR_Hispanic")

  combined_lrs<-left_join(combined_lrs,combined_lrs2)
  combined_lrs<-left_join(combined_lrs,combined_lrs3)
  combined_lrs<-left_join(combined_lrs,combined_lrs4)
  combined_lrs<-left_join(combined_lrs,combined_lrs5)
  combined_lrs<-left_join(combined_lrs,combined_lrs6)
  names(combined_lrs)[names(combined_lrs) == 'LR_relationship_tested'] <- 'L_R_relationship_tested'
  print(head(combined_lrs))

  combined_lrs<-combined_lrs %>% pivot_longer(
                             cols = starts_with("LR_"),
                             names_to = "population_tested",
                             names_prefix = "LR_",
                             values_to = "LR_population_tested",
                             values_drop_na = TRUE)

  names(combined_lrs)[names(combined_lrs) == 'L_R_relationship_tested'] <- 'LR_relationship_tested'
  View(combined_lrs)
  return(combined_lrs)
}

#######################################################################################################################

log_message("Processing individuals genotypes...")

processing_time <- system.time({
  processed_genotypes <- kinship_calculation(individuals_genotypes, kinship_matrix,df_allelefreq)
})
print(head(processed_genotypes))
log_message(paste("Processed individuals genotypes in", processing_time["elapsed"], "seconds."))
processed_genotypes <- as.data.table(processed_genotypes[, c("alleles_ind1","alleles_ind2","which_1st","A","B","C","D","marker"):=NULL])
print(head(processed_genotypes))
# Calculate combined likelihood ratios
log_message("Calculating combined likelihood ratios...")
combined_lrs_time <- system.time({
  combined_lrs <- calculate_combined_lrs(processed_genotypes, loci_lists)
})

rm(loci_lists)
log_message(paste("Calculated combined likelihood ratios in", combined_lrs_time["elapsed"], "seconds."))
#combined_lrs$LR_population_known<-as.numeric(combined_lrs$LR_population_known)

# Save results to CSV
log_message("Writing results to compressed files...")

setwd(output_dir)
write.table(combined_lrs,gzfile(paste0(output_file2,".gz")),append=FALSE)
setwd(t)
log_message("LR calculations completed.")
