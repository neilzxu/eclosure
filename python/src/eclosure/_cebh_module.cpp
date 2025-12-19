// Python bindings for the shared closed e-BH core in r/eclosure/src/closedebh.cpp.

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

#include <vector>

// Forward declarations of the C++ entry points we bind.
int closedeBH_cpp(const std::vector<double>& evalues,
                  double alpha = 0.05,
                  const std::vector<int>& set = {},
                  bool approximate = false);
int largest_mean_consistent_k(const std::vector<double>& evalues,
                              double alpha = 0.05,
                              bool approximate = false);
bool is_subset_mean_consistent(const std::vector<double>& evalues,
                               const std::vector<int>& subset,
                               double alpha = 0.05);

namespace py = pybind11;

PYBIND11_MODULE(_cebh, m) {
  m.doc() = "Pybind11 bindings for the closed e-BH core routines";

  m.def(
      "closede_bh_cpp",
      [](const std::vector<double>& evalues, double alpha, bool approximate) {
        return closedeBH_cpp(evalues, alpha, {}, approximate);
      },
      py::arg("evalues"),
      py::arg("alpha") = 0.05,
      py::arg("approximate") = false);

  m.def(
      "largest_mean_consistent_k",
      [](const std::vector<double>& evalues, double alpha, bool approximate) {
        return largest_mean_consistent_k(evalues, alpha, approximate);
      },
      py::arg("evalues"),
      py::arg("alpha") = 0.05,
      py::arg("approximate") = false);

  m.def(
      "is_subset_mean_consistent",
      [](const std::vector<double>& evalues,
         const std::vector<int>& subset,
         double alpha) {
        return is_subset_mean_consistent(evalues, subset, alpha);
      },
      py::arg("evalues"),
      py::arg("subset"),
      py::arg("alpha") = 0.05);
}
