closedBY = function(p, alpha) {
  m = length(p)
  p = sort(p)
  
  # Prepare all e-contributions for all |S| using BY-calibrated e-values
  es = t(vapply(
    seq_len(m),
    function(s) BY_calibrate_evalues(p, alpha, K = s),
    numeric(m)
  ))
  # pre-calculate cumulative sums
  cumes = t(apply(es, 1, cumsum))
  # print('new closedBY')
  # cumes which each row divided by it's row index
  cumes_norm = es / seq_len(m)
  # print(cumes_norm[1,])
  # print(p)
  # print(alpha / cumsum(1/1:m))
  
  # Check values or $r=|R|$ from r=m down
  r=m
  ready = FALSE
  while (!ready) {
    ready = TRUE
    # s=|S|
    for (s in m:1) {    
      # rs = |R cap S|
      # note that the size of rs is constrained by r and s
      for (rs in max(1,r+s-m):min(r,s)) {
        if (rs == r)
          ES = cumes[s,r] + cumes[s,m] - cumes[s,m- (s - r)]
        else
          ES = cumes[s,r] - cumes[s,r-rs] + cumes[s,m] - cumes[s,m-s+rs]

        critical = rs / (r * alpha)
        # stop this r as soon as we violate ES >= critical
        # prevent rounding errors
        if (ES / s < critical - 1e-10) {
          ready = FALSE
          r = r - 1
          break
        }
      }
      if (!ready) break
    }
    if (r==0) ready=TRUE
  }
  r
}

closedBY_old = function(p, alpha) {
  m = length(p)
  p = sort(p)
  
  # Prepare all e-contributions for all |S|
  longp = rep(p, each=m)
  longs = rep(1:m, m)
  hs = cumsum(1/1:m)
  longhs = rep(hs, m)
  #print('longp, longs, longhs:')
  #print(longp)
  #print(longs)
  #print(longhs)
  # put in a matrix
  #print('synthetic long p')
  #print(matrix(longp, m, m))
  #print(matrix(longhs, m, m))
  #print(matrix(longs, m, m))
  # columns are hypotheses, rows are 1:|S|
  es_vec = (longp*longhs <= alpha) / (alpha * ceiling(longs*longhs*longp/alpha))
  #print('es vec')
  #print(es_vec)
  es = matrix(es_vec, m, m)
  # pre-calculate cumulative sums
  #print('old closedBY')
  
  cumes = t(apply(es, 1, cumsum))
  #print(es[1,])
  
  # Check values or $r=|R|$ from r=m down
  r=m
  ready = FALSE
  while (!ready) {
    ready = TRUE
    # s=|S|
    for (s in m:1) {    
      # rs = |R cap S|
      # note that the size of rs is constrained by r and s
      for (rs in max(1,r+s-m):min(r,s)) {
        if (rs == r)
          ES = cumes[s,r] + cumes[s,m] - cumes[s,m-s+r]
        else
          ES = cumes[s,r] - cumes[s,r-rs] + cumes[s,m] - cumes[s,m-s+rs]
        critical = rs/(r*alpha) 
        # stop this r as soon as we violate ES >= critical
        # prevent rounding errors
        if (ES < critical - 1e-10) {
          ready = FALSE
          r = r - 1
          break
        }
      }
      if (!ready) break
    }
    if (r==0) ready=TRUE
  }
  r
}


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

