=======================
normal-inverse-gaussian
=======================

Description
===========

This library implements the cumulative distribution function of the normal inverse Gaussian (NIG) distribution.
The code is written in C++ and includes a Python package installable via pip.

C++ Build
=========

Requirements: a C++17 compiler and CMake >= 3.15.

.. code-block:: bash

   cmake -B build -DCMAKE_BUILD_TYPE=Release
   cmake --build build

This produces:

- ``build/libnig.a`` — static library to link against from C++ code.
- ``build/nig_demo`` — standalone executable.

The public header is ``code/include/nig.hpp``. Link against ``libnig`` and
include the header to call ``nig_cdf`` directly from C++.

Python Package
==============

Requirements: a C++17 compiler, CMake >= 3.15, and pybind11.

.. code-block:: bash

   pip install .

This builds and installs the ``nig`` package. No manual CMake invocation is
needed; scikit-build-core handles the compilation transparently.

Usage
=====

Two interfaces are provided.

**Wikipedia parametrisation** (α, β, μ, δ):

.. code-block:: python

   from nig import nig_cdf

   nig_cdf(x=2.0, alpha=2.0, beta=-0.4, mu=1.75, delta=2.0)

**SciPy-compatible parametrisation** — drop-in replacement for
``scipy.stats.norminvgauss``:

.. code-block:: python

   from nig import norminvgauss

   # unbound call (same signature as scipy.stats.norminvgauss.cdf)
   norminvgauss.cdf(x=2.0, a=4.0, b=-0.8, loc=1.75, scale=2.0)

   # frozen distribution
   dist = norminvgauss(a=4.0, b=-0.8, loc=1.75, scale=2.0)
   dist.cdf(2.0)

The relationship between the two parametrisations is
``a = alpha * delta``, ``b = beta * delta``, ``loc = mu``, ``scale = delta``.


Citation
========

If you use the library, please cite the paper https://arxiv.org/abs/2502.16015::

  @article{Navas-Palencia2025NIG,
    title     = {On the computation of the cumulative distribution function of the Normal
  Inverse Gaussian distribution},
    author    = {Navas-Palencia, G.},
    year      = {2025},
    eprint    = {2502.16015},
    archivePrefix = {arXiv},
    primaryClass = {math.NA},
    volume    = {abs/2502.16015},
    url       = {http://arxiv.org/abs/2502.16015},
  }