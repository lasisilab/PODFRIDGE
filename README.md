# PODFRIDGE

A [workflowr][] project.

[workflowr]: https://github.com/workflowr/workflowr

## Quick Start: Install All Required R Packages

Before running any analysis or building the site, you must install all required R packages. This project provides an automated script to do this for you.

**To install and update all required packages:**

### Option 1: In the R Console
```r
source("install_project_packages.R")
```

### Option 2: In the Terminal
```sh
Rscript install_project_packages.R
```

This script will:
- Scan all `.R` and `.Rmd` files for `library()` calls
- Install any missing packages
- Update all detected packages to the latest version
- Ensure the `workflowr` package is installed and up-to-date

**You should run this script before using workflowr or running any analysis.**
