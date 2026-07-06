#include <pybind11/pybind11.h>
#include "nig.hpp"

namespace py = pybind11;

PYBIND11_MODULE(_nig, m)
{
    m.doc() = "Normal Inverse Gaussian cumulative distribution function";

    m.def(
        "nig_cdf",
        &nig_cdf,
        py::arg("x"),
        py::arg("alpha"),
        py::arg("beta"),
        py::arg("mu"),
        py::arg("delta"),
        R"pbdoc(
            Compute the Normal Inverse Gaussian (NIG) CDF.

            Parameters
            ----------
            x : float
                Evaluation point.
            alpha : float
                Tail heaviness (alpha > |beta|).
            beta : float
                Asymmetry parameter.
            mu : float
                Location parameter.
            delta : float
                Scale parameter (delta > 0).

            Returns
            -------
            float
                CDF value at x.
        )pbdoc"
    );
}
