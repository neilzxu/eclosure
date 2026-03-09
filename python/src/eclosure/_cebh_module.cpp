// Python bindings for the guarded C++ routines in r/eclosure/src/ceBH.cpp and cBY.cpp.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <algorithm>
#include <numeric>
#include <stdexcept>
#include <vector>

// Forward declarations of the vector-based entry points compiled from ceBH.cpp.
int largestmeanconsistent_cpp(const std::vector<double>& e, double alpha);
int largestmeanconsistent_approximate_cpp(const std::vector<double>& e, double alpha);
int eBH_cpp(const std::vector<double>& e, double alpha);
bool meanconsistent_wrapper_cpp(const std::vector<double>& cum_e, int r, double alpha);

// Forward declarations of the vector-based entry points compiled from cBY.cpp.
struct cBYCheckResult {
  bool res;
  int at_s;
};

int BY_cpp(const std::vector<double>& p,
           const std::vector<double>& harmonic,
           double alpha);
cBYCheckResult cBY_check_cpp(const std::vector<double>& p,
                             int r,
                             const std::vector<double>& harmonic,
                             int warm_s,
                             double alpha);
int largestcBYsignificant_cpp(const std::vector<double>& p, double alpha);
int largestcBYsignificant_approximate_cpp(const std::vector<double>& p,
                                          double alpha);

namespace py = pybind11;

namespace {

std::vector<double> build_cum_e_for_subset(const std::vector<double>& evalues,
                                           const std::vector<int>& subset) {
  std::vector<bool> in_subset(evalues.size(), false);
  for (int idx : subset) {
    if (idx < 0 || idx >= static_cast<int>(evalues.size())) {
      throw std::out_of_range("subset index out of range");
    }
    in_subset[idx] = true;
  }

  std::vector<double> e_in;
  std::vector<double> e_out;
  e_in.reserve(subset.size());
  e_out.reserve(evalues.size() - subset.size());

  for (int i = 0; i < static_cast<int>(evalues.size()); ++i) {
    if (in_subset[i]) {
      e_in.push_back(evalues[i]);
    } else {
      e_out.push_back(evalues[i]);
    }
  }

  std::sort(e_in.begin(), e_in.end());
  std::sort(e_out.begin(), e_out.end());

  std::vector<double> cum_e(evalues.size() + 1, 0.0);
  int pos = 1;
  for (double value : e_out) {
    cum_e[pos] = cum_e[pos - 1] + value;
    ++pos;
  }
  for (double value : e_in) {
    cum_e[pos] = cum_e[pos - 1] + value;
    ++pos;
  }
  return cum_e;
}

}  // namespace

PYBIND11_MODULE(_cebh, m) {
  m.doc() = "Pybind11 bindings for the eClosure C++ routines";

  m.def(
      "closede_bh_cpp",
      [](const std::vector<double>& evalues, double alpha, bool approximate) {
        if (approximate) {
          return largestmeanconsistent_approximate_cpp(evalues, alpha);
        }
        return largestmeanconsistent_cpp(evalues, alpha);
      },
      py::arg("evalues"),
      py::arg("alpha") = 0.05,
      py::arg("approximate") = false);

  m.def(
      "largest_mean_consistent_k",
      [](const std::vector<double>& evalues, double alpha, bool approximate) {
        if (approximate) {
          return largestmeanconsistent_approximate_cpp(evalues, alpha);
        }
        return largestmeanconsistent_cpp(evalues, alpha);
      },
      py::arg("evalues"),
      py::arg("alpha") = 0.05,
      py::arg("approximate") = false);

  m.def(
      "is_subset_mean_consistent",
      [](const std::vector<double>& evalues,
         const std::vector<int>& subset,
         double alpha) {
        if (subset.empty()) {
          return true;
        }
        std::vector<double> cum_e = build_cum_e_for_subset(evalues, subset);
        return meanconsistent_wrapper_cpp(
            cum_e, static_cast<int>(subset.size()), alpha);
      },
      py::arg("evalues"),
      py::arg("subset"),
      py::arg("alpha") = 0.05);

  // Raw routines exposed for the upcoming API migration.
  m.def("largestmeanconsistent_cpp",
        &largestmeanconsistent_cpp,
        py::arg("evalues"),
        py::arg("alpha"));
  m.def("largestmeanconsistent_approximate_cpp",
        &largestmeanconsistent_approximate_cpp,
        py::arg("evalues"),
        py::arg("alpha"));
  m.def("eBH_cpp", &eBH_cpp, py::arg("evalues"), py::arg("alpha"));
  m.def("meanconsistent_wrapper_cpp",
        &meanconsistent_wrapper_cpp,
        py::arg("cum_e"),
        py::arg("r"),
        py::arg("alpha"));
  m.def("BY_cpp", &BY_cpp, py::arg("pvalues"), py::arg("harmonic"), py::arg("alpha"));
  m.def("largestcBYsignificant_cpp",
        &largestcBYsignificant_cpp,
        py::arg("pvalues"),
        py::arg("alpha"));
  m.def("largestcBYsignificant_approximate_cpp",
        &largestcBYsignificant_approximate_cpp,
        py::arg("pvalues"),
        py::arg("alpha"));
  m.def(
      "cBY_check_cpp",
      [](const std::vector<double>& pvalues,
         int r,
         const std::vector<double>& harmonic,
         int warm_s,
         double alpha) {
        cBYCheckResult result = cBY_check_cpp(pvalues, r, harmonic, warm_s, alpha);
        py::dict out;
        out["res"] = result.res;
        out["at_s"] = result.at_s;
        return out;
      },
      py::arg("pvalues"),
      py::arg("r"),
      py::arg("harmonic"),
      py::arg("warm_s"),
      py::arg("alpha"));
}
