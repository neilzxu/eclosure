#ifndef E_CLOSURE_NO_RCPP
#include <Rcpp.h>
using namespace Rcpp;
#endif

#include <algorithm>
#include <cmath>
#include <limits>
#include <queue>
#include <stdexcept>
#include <vector>
#include <stdexcept>

using std::vector;

/**********************************************************************
 * meanconsistent  (the updated version with early termination)
 **********************************************************************/
bool meanconsistent(
    const vector<double>& f,
    const vector<double>& g,
    double alpha,
    int r,
    int m
)
{
  const double tol = 1e-12;
  if (r == 0) return true;
  const int Jmax = m - r;

  auto h = [&](int i, int j) {
    double fi = f[i];
    double gj = g[j];
    double coeff = static_cast<double>(i)/(r*alpha);
    return fi + gj - (i + j)*coeff;
  };

  // Start with the "root" of the implicit balanced binary tree.
  // The original code pushed r first, which corresponds to the root interval [r, r].
  // Then it pushed [1, r-1]. We emulate that:
  std::queue<std::pair<int,int>> qi;
  qi.push({r, r});
  if (r > 1) qi.push({1, r - 1});

  while (!qi.empty()) {
    auto [L, R] = qi.front();
    qi.pop();
    if (L > R) continue;

    int mid = (L + R) / 2;
    int i = mid;

    // ---- convex search in j ----
    int left = 0;
    int right = Jmax;

    while (left + 3 <= right) {
      int m1 = left + (right - left) / 3;
      int m2 = right - (right - left) / 3;

      double h1 = h(i, m1);
      if (h1 < -tol) return false;

      double h2 = h(i, m2);
      if (h2 < -tol) return false;

      if (h1 < h2)
        right = m2 - 1;
      else
        left = m1 + 1;
    }

    for (int j = left; j <= right; ++j)
      if (h(i, j) < -tol) return false;

    // ---- if we get here, mid passes; push its children ----
    if (L <= mid - 1) qi.push({L, mid - 1});
    if (mid + 1 <= R) qi.push({mid + 1, R});
  }

  return true; // all (i,j) checked and nonnegative
}

/**********************************************************************
 * Helper: cumulative sums
 **********************************************************************/
static void build_cumsum_set(
    const vector<double>& e_sorted,
    const vector<int>& idx_set,
    vector<double>& f
) {
  int r = idx_set.size();
  f.resize(r+1);
  f[0] = 0.0;
  for (int i = 0; i < r; i++)
    f[i+1] = f[i] + e_sorted[idx_set[i]];
}

static void build_cumsum_complement(
    const vector<double>& e_sorted,
    const vector<int>& idx_set,
    vector<double>& g
) {
  int m = e_sorted.size();
  int r = idx_set.size();
  g.resize(m-r+1);
  g[0] = 0.0;

  vector<int> comp;
  comp.reserve(m-r);
  int p = 0;
  for (int i = 0; i < m; i++) {
    if (p < r && idx_set[p] == i) p++;
    else comp.push_back(i);
  }
  for (int j = 0; j < (int)comp.size(); j++)
    g[j+1] = g[j] + e_sorted[comp[j]];
}

/**********************************************************************
 * Core helpers reusable across R and Python builds
 **********************************************************************/
bool is_subset_mean_consistent(const vector<double>& evalues,
                               const vector<int>& subset,
                               double alpha) {
  int m = static_cast<int>(evalues.size());
  if (m == 0 || subset.empty()) return true;

  vector<double> e_sorted(evalues.begin(), evalues.end());
  std::sort(e_sorted.begin(), e_sorted.end());

  vector<std::pair<double, int>> val_idx(m);
  for (int i = 0; i < m; i++) val_idx[i] = {evalues[i], i};
  std::sort(val_idx.begin(), val_idx.end(),
            [](auto &a, auto &b){ return a.first < b.first; });

  vector<int> pos(m);
  for (int i = 0; i < m; i++) pos[val_idx[i].second] = i;

  vector<int> idx_set;
  idx_set.reserve(subset.size());
  for (auto xi : subset) {
    if (xi < 0 || xi >= m) {
      throw std::out_of_range("subset index out of range");
    }
    idx_set.push_back(pos[xi]);
  }
  std::sort(idx_set.begin(), idx_set.end());

  vector<double> f, g;
  build_cumsum_set(e_sorted, idx_set, f);
  build_cumsum_complement(e_sorted, idx_set, g);

  return meanconsistent(f, g, alpha, static_cast<int>(idx_set.size()), m);
}

int largest_mean_consistent_k(const vector<double>& evalues,
                              double alpha,
                              bool approximate) {
  int m = static_cast<int>(evalues.size());
  if (m == 0) return 0;

  vector<double> e_sorted(evalues.begin(), evalues.end());
  std::sort(e_sorted.begin(), e_sorted.end());

  vector<double> prefix(m+1);
  prefix[0] = 0.0;
  for (int i = 0; i < m; i++) prefix[i+1] = prefix[i] + e_sorted[i];

  auto check_k = [&](int k){
    int r = k;
    int s = m - k;
    vector<double> f(k+1), g(s+1);
    f[0] = 0.0;
    for (int i=1;i<=k;i++) f[i] = f[i-1] + e_sorted[m-k + (i-1)];
    for (int j=0;j<=s;j++) g[j] = prefix[j];
    return meanconsistent(f,g,alpha,r,m);
  };

  if (!approximate) {
    for (int k=m;k>=0;k--) if (check_k(k)) return k;
    return 0;
  }

  int lo=0, hi=m, best=0;
  while(lo<=hi){
    int mid = (lo+hi)/2;
    if(check_k(mid)){ best=mid; lo=mid+1; }
    else hi=mid-1;
  }

  return best;
}

/**********************************************************************
 * Main function exposed to R (and reused in Python via pybind)
 **********************************************************************/
#ifndef E_CLOSURE_NO_RCPP
// [[Rcpp::export]]
int closedeBH_cpp(NumericVector e,
                  double alpha = 0.05,
                  Nullable<IntegerVector> set = R_NilValue,
                  bool approximate = false)
#else
int closedeBH_cpp(const vector<double>& e,
                  double alpha = 0.05,
                  const vector<int>& set = {},
                  bool approximate = false)
#endif
{
  vector<double> evalues(e.begin(), e.end());
  int m = static_cast<int>(evalues.size());
  if (m == 0) return 0;

  // CASE 1: set is provided
#ifndef E_CLOSURE_NO_RCPP
  bool have_set = set.isNotNull();
#else
  bool have_set = !set.empty();
#endif
  if (have_set) {
#ifndef E_CLOSURE_NO_RCPP
    IntegerVector setR(set);
    vector<int> raw_set(setR.begin(), setR.end());
#else
    vector<int> raw_set(set.begin(), set.end());
#endif

    // Handle negative indexing like in R
    vector<int> final_set;
    if (std::any_of(raw_set.begin(), raw_set.end(), [](int x){ return x < 0; })) {
      for (int i = 0; i < m; i++)
        if (std::find(raw_set.begin(), raw_set.end(), -(i+1)) == raw_set.end())
          final_set.push_back(i);
    } else {
#ifndef E_CLOSURE_NO_RCPP
      for (auto x : raw_set) final_set.push_back(x-1); // 1-based → 0-based
#else
      for (auto x : raw_set) final_set.push_back(x);   // assume 0-based for non-R builds
#endif
    }

    bool ok = is_subset_mean_consistent(evalues, final_set, alpha);
    return ok ? 1 : 0;
  }

  // CASE 2: set missing → approximate or exact search
  return largest_mean_consistent_k(evalues, alpha, approximate);
}
