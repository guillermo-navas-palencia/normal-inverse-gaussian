import numpy.ctypeslib as npct
import os.path
import platform

from ctypes import c_double


# load library
libabspath = os.path.dirname(os.path.abspath(__file__))

lib = npct.load_library('_nig.so', libabspath)


lib.cpp_nig_cdf.restype = c_double
lib.cpp_nig_cdf.argtypes = (c_double, c_double, c_double, c_double, c_double)

def nig_cdf(
    x: float,
    alpha: float,
    beta: float,
    mu: float,
    delta: float
) -> float:

    return lib.cpp_nig_cdf(x, alpha, beta, mu, delta)