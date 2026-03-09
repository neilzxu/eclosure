harmonic_number = function(K) {
  if (K <= 0) {
    return(0)
  }
  sum(1 / seq_len(K))
}


BY_calibrate_evalues = function(pvalues, alpha, K = length(pvalues)) {
  ell_K = harmonic_number(K)
  evalues = rep(0, length(pvalues))
  incr = alpha / (K * ell_K)
  if (incr <= 0) {
    return(evalues)
  }
  mask = pvalues <= (alpha / ell_K)
  if (any(mask)) {
    denom = alpha * pmax(ceiling(pvalues[mask] / incr), 1)
    evalues[mask] = K / denom
  }
  evalues
}


ebh_discoveries = function(evalues, alpha) {
  m = length(evalues)
  if (m == 0) {
    return(0)
  }
  sorted = sort(evalues, decreasing = TRUE)
  thresholds = m / (alpha * seq_len(m))
  eligible = which(sorted >= thresholds)
  if (length(eligible) == 0) {
    return(0)
  }
  max(eligible)
}


closed_ebh_from_evalues = function(evalues, alpha) {
  m = length(evalues)
  if (m == 0) {
    return(0)
  }

  order_idx = order(evalues, decreasing = TRUE)
  sorted = evalues[order_idx]
  ebh_r = ebh_discoveries(sorted, alpha)

  tail_sum = function(k, r) {
    if (r <= 0) {
      return(0)
    }
    start = k - (r - 1)
    if (start < 1) {
      stop("Invalid range in tail_sum")
    }
    sum(sorted[start:k])
  }

  if (ebh_r == m) {
    return(m)
  }

  for (k in seq(m, ebh_r + 1, by = -1)) {
    valid = TRUE
    for (r_false in seq_len(k)) {
      for (m_delta in 0:(m - k)) {
        m_null = r_false + m_delta
        rejected_sum = tail_sum(k, r_false)
        nonreject_sum = tail_sum(m, m_delta)
        evalue = (rejected_sum + nonreject_sum) / m_null
        fdp = r_false / k
        if (evalue < (fdp / alpha) - 1e-10) {
          valid = FALSE
          break
        }
      }
      if (!valid) {
        break
      }
    }
    if (valid) {
      return(k)
    }
  }

  ebh_r
}


closed_eBH_cal = function(p, alpha) {
  m = length(p)
  if (m == 0) {
    return(0)
  }
  evalues = BY_calibrate_evalues(p, alpha, m)
  closed_ebh_from_evalues(evalues, alpha)
}
