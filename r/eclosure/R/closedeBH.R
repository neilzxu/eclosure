#' Closed e-BH Mean-Consistency Check and Set Size Finder
#'
#' @description
#' Implements the closed e-Benjamini–Hochberg (eBH) mean-consistency
#' condition for a given set of hypotheses or searches for the largest
#' mean-consistent set among the \eqn{k} largest e-values.
#'
#' The function either:
#' \itemize{
#'   \item checks whether a user-specified set is *mean-consistent*, or
#'   \item (if no set is supplied) finds a set size \eqn{k} such that the
#'         \eqn{k} largest e-values form a mean-consistent set.
#' }
#'
#' The underlying consistency check is performed by a fast C++ routine
#' based on convex minimization over the cumulative sums of e-values.
#'
#'
#' @details
#' Let \eqn{e_1,\dots, e_m} be e-values and let \eqn{\alpha \in (0,1]}.
#'
#' A set \eqn{S} of size \eqn{r} is called *mean-consistent* if, after
#' ordering the e-values in \(S\) and in its complement in ascending order,
#' the function
#'
#' \deqn{
#' h(i,j) = f_i + g_j - (i + j)\,\frac{i}{r\alpha}
#' }
#'
#' is non-negative for all \eqn{0 \le i \le r} and
#' \eqn{0 \le j \le m-r}. Here
#' \itemize{
#'   \item \eqn{f_i} is the cumulative sum of the smallest \eqn{i} e-values
#'         inside the set \eqn{S},
#'   \item \eqn{g_j} is the cumulative sum of the smallest \eqn{j} e-values
#'         outside the set \eqn{S}.
#' }
#'
#' The C++ routine guarantees early termination whenever a violating
#' \eqn{(i,j)} pair is found. For each fixed \eqn{i}, the function
#' \eqn{h(i,j)} is convex in \eqn{j}, and a discrete convex binary search
#' is used to find its minimum.
#'
#' ## Behaviour depending on `set`:
#'
#' - **If `set` is provided:**
#'   The function checks the mean-consistency of exactly this set and returns
#'   a logical value.
#'   The `approximate` argument is ignored.
#'
#' - **If `set` is not provided:**
#'   The function considers sets consisting of the \eqn{k} largest e-values
#'   for \eqn{k=0,\dots,m}.
#'
#'   - If `approximate = FALSE`, a linear search from \eqn{k=m} downward
#'     is used to find the largest mean-consistent set.
#'   - If `approximate = TRUE`, a bisection search is used.
#'     The returned \eqn{k} always yields a mean-consistent set but may not be
#'     the largest.
#'
#' The computational complexity is approximately:
#' \itemize{
#'   \item \eqn{O(m \log m)} when `set` is provided,
#'   \item \eqn{O(m^2 \log m)} when `set` is missing and `approximate = FALSE`,
#'   \item \eqn{O(m \log^2 m)} when `set` is missing and `approximate = TRUE`.
#' }
#'
#'
#' @param e A numeric vector of length \eqn{m} containing e-values.
#' @param alpha A number in \eqn{[0,1]}, the target FDR level.
#'   Default is \code{0.05}.
#' @param set Optional integer or logical vector specifying a subset of
#'   \code{1:m} to check for mean-consistency. If provided, the function
#'   returns a logical value. If \code{NULL} (default), the function searches
#'   for a mean-consistent set of largest possible size.
#' @param approximate Logical (default \code{FALSE}). Ignored if `set` is
#'   provided.
#'   If `TRUE`, bisection search is used to find a (not necessarily maximal)
#'   mean-consistent set.
#'   If `FALSE`, a descending linear search from \code{m} is used to
#'   find the *largest* mean-consistent set.
#' @param e Numeric vector of e-values.
#' @param alpha Numeric, significance level (default 0.05).
#' @param set Optional set (logical, positive integer, or negative integer).
#' @param approximate Logical, whether to use approximate bisection (default TRUE).
#' @return If `set` is provided: TRUE/FALSE. Otherwise: size of the largest mean-consistent set.
#' @export
#'
closedeBH <- function(e, alpha=0.05, set=NULL, approximate=FALSE) {
  if (!is.null(set)) {
    if (is.logical(set)) {
      if (length(set) != length(e))
        stop("Logical 'set' must be same length as e")
      set <- which(set)
    }
  }
  res <- closedeBH_cpp(e, alpha, set, approximate)
  if (!is.null(set)) return(as.logical(res))
  return(res)
}
#'
#' @return
#' - If `set` is provided: a single logical value, indicating whether
#'   the set is mean-consistent.
#' - If `set` is not provided: a single integer giving the size of a
#'   mean-consistent set (largest such set if `approximate = FALSE`).
#'
#'
#' @examples
#' set.seed(1)
#' e <- rexp(100)
#'
#' # Check a given set
#' closedeBH(e, set = 1:10)
#'
#' # Search for a consistent set size (approximate bisection)
#' closedeBH(e)
#'
#' # Exact search for the largest mean-consistent set
#' closedeBH(e, approximate = FALSE)
#'
#' @useDynLib eClosure, .registration = TRUE
#' @importFrom Rcpp sourceCpp
