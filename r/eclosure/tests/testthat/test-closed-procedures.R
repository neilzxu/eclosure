test_that("closed eBH dominates eBH on packaged real e-values", {
  data_dir <- system.file("extdata", "evalues", package = "eClosure")
  skip_if(data_dir == "")
  csvs <- list.files(data_dir, pattern = "[.]csv$", full.names = TRUE)
  skip_if(length(csvs) == 0)

  for (alpha in c(0.05, 0.1)) {
    for (path in csvs) {
      print(path)
      df <- utils::read.csv(path)
      expect_true("evalue" %in% names(df))
      evalues <- na.omit(df$evalue)
      closed_k <- closedeBH(evalues, alpha = alpha, approximate = FALSE)
      ebh_k <- ebh_discoveries(evalues, alpha = alpha)
      
      expect_gte(
        closed_k,
        ebh_k,
        label = sprintf("File %s at alpha %.2f", basename(path), alpha)
      )
    }
  }
})

test_that("closed eBH dominates eBH on simulated e-values", {
  set.seed(123)
  evalues <- rexp(128, rate = 1 / 1.5)
  for (alpha in c(0.05, 0.1)) {
    closed_k <- closedeBH(evalues, alpha = alpha, approximate = FALSE)
    ebh_k <- ebh_discoveries(evalues, alpha = alpha)
    expect_gte(closed_k, ebh_k)
  }
})

test_that("closed BY discoveries include BY on real p-values", {
  data_dir <- system.file("extdata", "pvalues", package = "eClosure")
  skip_if(data_dir == "")
  csvs <- list.files(data_dir, pattern = "[.]csv$", full.names = TRUE)
  skip_if(length(csvs) == 0)

  for (alpha in c(0.05, 0.1)) {
    for (path in csvs) {
      df <- utils::read.csv(path)
      expect_true("pvalue" %in% names(df))
      pv <- df$pvalue
      closed_by <- closedBY(pv, alpha)
      by_count <- sum(stats::p.adjust(pv, method = "BY") <= alpha)
      expect_gte(
        closed_by,
        by_count,
        label = sprintf("File %s at alpha %.2f", basename(path), alpha)
      )
    }
  }
})

test_that("closed BY dominates BY on simulated p-values", {
  set.seed(321)
  pvals <- sort(runif(200))
  for (alpha in c(0.05, 0.1)) {
    closed_by <- closedBY(pvals, alpha)
    by_count <- sum(stats::p.adjust(pvals, method = "BY") <= alpha)
    expect_gte(closed_by, by_count)
  }
})
