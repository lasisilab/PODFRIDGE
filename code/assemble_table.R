library(data.table)

args <- commandArgs(trailingOnly = TRUE)
str(args)
cat(args, sep = "\n")
job_name<-args[1]

# Helper function for logging
log_message <- function(message) {
  cat(paste0("[", Sys.time(), "] ", message, "\n"))
}
log_message("Start table assembly")

t<-getwd()
print(job_name)
print(paste0(t,"/data/sims"))
setwd(paste0(t,"/data/sims"))
tabs<-list.files(pattern = "processed_genotypes.csv", recursive = TRUE)
final_results<-do.call(rbind, lapply(tabs, read.csv))

array_task<-rep(1:length(tabs),each=(nrow(final_results)/length(tabs)))

output_dir <- file.path(t,"data", paste0("simulation_", job_name))
dir.create(output_dir, recursive = TRUE)

output_file <- file.path(output_dir, "processed_genotypes.csv")
fwrite(final_results, output_file,append=TRUE)
log_message("Table assembled.")
