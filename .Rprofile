## This makes sure that R loads the workflowr package
## automatically, everytime the project is loaded in RStudio
if (Sys.getenv("RSTUDIO") == "1") {
  if (requireNamespace("workflowr", quietly = TRUE)) {
    message("Loading .Rprofile for the current workflowr project")
    library("workflowr")
  } else {
    message("workflowr package not installed, please run install.packages(\"workflowr\") to use the workflowr functions")
  }
}
