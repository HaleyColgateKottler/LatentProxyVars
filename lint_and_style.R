library(lintr)
library(styler)

# Function to lint R scripts
lint_dir <- function(dir) {
  files <- list.files(dir, pattern = "\\.R$", full.names = TRUE)
  lapply(files, lintr::lint)
}

# Function to style R scripts
style_dir <- function(dir) {
  files <- list.files(dir, pattern = "\\.R$", full.names = TRUE)
  lapply(files, styler::style_file)
}

# Define the project directory (replace with your project's path)
project_dir <- getwd()

# Run lintr and styler on the project directory
lint_dir(project_dir)
style_dir(project_dir)
style_dir('tests')
lint("main.R")
