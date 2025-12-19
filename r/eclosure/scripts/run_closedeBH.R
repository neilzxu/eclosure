#!/usr/bin/env Rscript

suppressPackageStartupMessages({
  if (!requireNamespace("pkgload", quietly = TRUE)) {
    stop("The pkgload package is required to run this helper script.")
  }
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 5) {
  stop("Usage: run_closedeBH.R <pkg_dir> <e_path> <alpha> <approximate> <set_path_or_NULL>")
}

pkg_dir <- normalizePath(args[1], mustWork = TRUE)
e_path <- normalizePath(args[2], mustWork = TRUE)
alpha <- as.numeric(args[3])
approximate <- as.logical(args[4])
set_path <- args[5]

pkgload::load_all(pkg_dir, quiet = TRUE)

read_vector <- function(path) {
  if (!file.exists(path) || file.info(path)$size == 0) {
    return(numeric(0))
  }
  scan(path, quiet = TRUE)
}

evalues <- read_vector(e_path)

set_values <- NULL
if (!identical(set_path, "NULL")) {
  set_values <- as.integer(read_vector(normalizePath(set_path, mustWork = TRUE)))
}

result <- eClosure::closedeBH(evalues,
                              alpha = alpha,
                              set = set_values,
                              approximate = approximate)

if (is.logical(result)) {
  cat(ifelse(isTRUE(result), "TRUE", "FALSE"))
} else {
  cat(result)
}
