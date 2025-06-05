# install_project_packages.R
# This script scans all .Rmd and .R files in the project for required R packages and installs any that are missing.

get_packages_from_file <- function(filepath) {
  lines <- readLines(filepath, warn = FALSE)
  pkgs <- c()
  # Only find library() calls (not require or ::)
  lib_matches <- regmatches(lines, gregexpr("library\\((['\"]?)([A-Za-z0-9_.]+)\\1\\)", lines, perl=TRUE))
  for (m in lib_matches) {
    for (call in m) {
      # Extract just the package name (second capture group)
      pkg <- sub("library\\((['\"]?)([A-Za-z0-9_.]+)\\1\\)", "\\2", call, perl=TRUE)
      pkgs <- c(pkgs, pkg)
    }
  }
  unique(pkgs)
}

# Get all R and Rmd files, but only process regular files (not directories)
get_r_files <- function(path = ".") {
  files <- list.files(path, pattern = "\\.(R|r|Rmd|rmd)$", recursive = TRUE, full.names = TRUE)
  files[file.info(files)$isdir == FALSE]
}

r_files <- get_r_files()

all_pkgs <- unique(unlist(lapply(r_files, get_packages_from_file)))

# Remove base R and recommended packages (optional, can be customized)
base_pkgs <- rownames(installed.packages(priority = "base"))
recommended_pkgs <- rownames(installed.packages(priority = "recommended"))
all_pkgs <- setdiff(all_pkgs, c(base_pkgs, recommended_pkgs, ""))

# Check which packages are not installed
not_installed <- setdiff(all_pkgs, rownames(installed.packages()))

# Set a default CRAN mirror for non-interactive use
options(repos = c(CRAN = "https://cloud.r-project.org"))

if (length(not_installed) > 0) {
  cat("Installing missing packages:\n")
  print(not_installed)
  install.packages(not_installed)
} else {
  cat("All required packages are already installed.\n")
}

# Always ensure workflowr is installed and up-to-date
if (!requireNamespace("workflowr", quietly = TRUE)) {
  install.packages("workflowr")
} else {
  # Optionally update workflowr if not the latest version
  tryCatch({
    update.packages(oldPkgs = "workflowr", ask = FALSE)
  }, error = function(e) {
    message("Could not update workflowr: ", e$message)
  })
}

# Optionally update all detected packages if already installed
if (length(all_pkgs) > 0) {
  update.packages(oldPkgs = all_pkgs, ask = FALSE)
}

cat("\nAll detected required packages:\n")
print(sort(all_pkgs))
