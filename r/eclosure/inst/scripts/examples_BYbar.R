rm(list = ls())

alpha_levels <- c(0.05, 0.1)

source("closedBY.R")

suppressPackageStartupMessages({
  library(cherry)
  library(fdrtool)
})

# ------------------------------
# Dataset loaders
# ------------------------------

load_apsac <- function() {
  c(
    0.0001, 0.0004, 0.0019, 0.0095, 0.0201,
    0.0278, 0.0298, 0.0344, 0.0459, 0.3240,
    0.4262, 0.5719, 0.6528, 0.7590, 1.0000
  )
}

load_naep <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data_env <- new.env()
      data("NAEP", package = "cherry", envir = data_env)
      cache <<- as.numeric(data_env$NAEP) + 10^(-8)
    }
    cache
  }
})

load_padjust <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      old_seed_exists <- exists(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
      if (old_seed_exists) {
        old_seed <- get(".Random.seed", envir = .GlobalEnv, inherits = FALSE)
        on.exit(assign(".Random.seed", old_seed, envir = .GlobalEnv), add = TRUE)
      } else {
        on.exit(rm(".Random.seed", envir = .GlobalEnv), add = TRUE)
      }
      set.seed(123)
      x <- rnorm(50, mean = c(rep(0, 25), rep(3, 25)))
      cache <<- 2 * pnorm(sort(-abs(x)))
    }
    cache
  }
})

load_pvalues <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data_env <- new.env()
      data("pvalues", package = "fdrtool", envir = data_env)
      cache <<- as.numeric(data_env$pvalues)
    }
    cache
  }
})

load_vandevijver <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      data_env <- new.env()
      load("r/eclosure/inst/extdata/pvalues_tutorial.RData", envir = data_env)
      cache <<- as.numeric(data_env$ps)
    }
    cache
  }
})

load_golub <- local({
  cache <- NULL
  function() {
    if (is.null(cache)) {
      leukemia_big <- read.csv(
        "http://hastie.su.domains/CASI_files/DATA/leukemia_big.csv",
        check.names = FALSE
      )
      ix <- c(1:20, 35:61)
      cache <<- apply(
        leukemia_big,
        1,
        function(row) t.test(row[ix], row[-ix], equal.var = TRUE)$p.value
      )
    }
    cache
  }
})

datasets <- list(
  list(
    name = "APSAC",
    loader = load_apsac,
    reference = "\\cite{benjamini_controlling_false_1995}"
  ),
  list(
    name = "NAEP",
    loader = load_naep,
    reference = "\\cite{benjamini2000adaptive}"
  ),
  list(
    name = "PADJUST",
    loader = load_padjust,
    reference = "Example in \\texttt{p.adjust} R function"
  ),
  list(
    name = "PVALUES",
    loader = load_pvalues,
    reference = "Data in \\texttt{fdrtool} R package"
  ),
  list(
    name = "VANDEVIJVER",
    loader = load_vandevijver,
    reference = "\\cite{goeman2014multiple}"
  ),
  list(
    name = "GOLUB",
    loader = load_golub,
    reference = "\\citet{efron_computer_age_2021}"
  )
)




compute_method_counts <- function(pvalues, alpha) {
  list(
    BH = sum(p.adjust(pvalues, method = "BH") <= alpha),
    BY = sum(p.adjust(pvalues, method = "BY") <= alpha),
    #cBY_old = closedBY_old(pvalues, alpha),
    cBY = closedBY(pvalues, alpha),
    closed_eBH_cal = closed_eBH_cal(pvalues, alpha)
  )
}


results <- lapply(datasets, function(ds) {
  pvalues <- ds$loader()
  m <- length(pvalues)
  row <- data.frame(
    Dataset = ds$name,
    Hypotheses = m,
    stringsAsFactors = FALSE
  )

  for (alpha in alpha_levels) {
    key <- sprintf("%.2f", alpha)
    counts <- compute_method_counts(pvalues, alpha)
    row[[paste0("BH_", key)]] <- counts$BH
    row[[paste0("BY_", key)]] <- counts$BY
    #row[[paste0("cBYold_", key)]] <- counts$cBY_old
    row[[paste0("cBY_", key)]] <- counts$cBY
    row[[paste0("cEBH_", key)]] <- counts$closed_eBH_cal
  }

  row$Reference <- ds$reference
  row
})

results_df <- do.call(rbind, results)

print(results_df)


# ------------------------------
# LaTeX table
# ------------------------------

format_int <- function(x) {
  formatC(as.integer(round(x)), format = "d", big.mark = ",")
}

bold_if <- function(value, condition) {
  if (condition) {
    sprintf("\\textbf{%s}", value)
  } else {
    value
  }
}

alpha_keys <- sprintf("%.2f", alpha_levels)

cat("\\begin{table}[htbp]\n")
cat("\\centering\n")
cat("\\caption{Number of discoveries by BH, the \\citet{benjamini_false_discovery_2005} procedure (BY), our new procedure $\\cBY$ (closed BY), and closed eBH + calibrator across a variety of real datasets for FDR thresholds $\\alpha$. Settings where $\\cBY$ makes more discoveries than BY are highlighted in bold.}\n")
cat("\\label{tab:by_real_data_discoveries}\n")
cat("{\\setlength{\\tabcolsep}{3pt}\n")
cat("\\begin{tabular}{lccccccccccccc}\n")
cat("\\toprule\n")
cat("Dataset & \\makecell{\\# of\\\\hypotheses} & \\multicolumn{5}{c}{$\\alpha = 0.05$} & \\multicolumn{5}{c}{$\\alpha = 0.1$} & Reference \\\\n")
cat("\\cmidrule(lr){3-7} \\cmidrule(lr){8-12}\n")
cat(" &  & BH & BY & $\\cBY$ (old) & $\\cBY$ (ours) & closed eBH + cal & BH & BY & $\\cBY$ (old) & $\\cBY$ (ours) & closed eBH + cal &  \\\\\n")
cat("\\midrule\n")

for (i in seq_len(nrow(results_df))) {
  row <- results_df[i, , drop = FALSE]

  cells <- c(
    row[["Dataset"]],
    format_int(row[["Hypotheses"]])
  )

  for (key in alpha_keys) {
    bh_val <- row[[paste0("BH_", key)]]
    by_val <- row[[paste0("BY_", key)]]
    #cby_old_val <- row[[paste0("cBYold_", key)]]
    cby_val <- row[[paste0("cBY_", key)]]
    cebh_val <- row[[paste0("cEBH_", key)]]

    cells <- c(
      cells,
      format_int(bh_val),
      format_int(by_val),
      #format_int(cby_old_val),
      bold_if(format_int(cby_val), cby_val > by_val),
      format_int(cebh_val)
    )
  }

  cells <- c(cells, row[["Reference"]])
  cat(paste(cells, collapse = " & "), " \\\\n", sep = "")
}

cat("\\bottomrule\n")
cat("\\end{tabular}\n")
cat("}\n")
cat("\\end{table}\n")
