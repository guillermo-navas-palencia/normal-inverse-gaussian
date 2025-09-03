import os.path
import time

from ctypes import c_double

import numpy as np
import numpy.ctypeslib as npct
import numpy.typing as npt
import pandas as pd


# load library
libabspath = os.path.dirname(os.path.abspath(__file__))

lib = npct.load_library('../_nig.so', libabspath)

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


def analytical_european_option_call(
    St: float,
    K: float,
    r: float,
    q: float,
    tau: float,
    alpha: float,
    beta: float,
    delta: float
) -> float:
    
    gamma = np.sqrt(alpha**2 - beta**2)
    omega = delta * (np.sqrt(alpha**2 - (beta + 1)**2) - gamma)
    k = np.log(St / K) + (r - q + omega) * tau

    cdf1 = nig_cdf(x=k, alpha=alpha, beta=-beta, mu=0.0, delta=delta * tau)
    cdf2 = nig_cdf(x=k, alpha=alpha, beta=-beta-1, mu=0.0, delta=delta * tau)

    return St * np.exp(-q * tau) * cdf2 - K * np.exp(-r * tau) * cdf1


if __name__ == '__main__':
    St = 3000
   
    r = 0.01
    q = 0.0

    alpha = 40.0
    beta = 0.0
    delta = 25.0

    time_init = time.perf_counter()

    d_calls = {}
    for tau in (1/12, 1/6, 1, 2):
        calls = []
        for K in range(2000, 5000):
            call = analytical_european_option_call(
                St=St, K=K, r=r, q=q, tau=tau, alpha=alpha, beta=beta, delta=delta
            )
            calls.append(call)
        d_calls[str(tau)] = calls

    time_elapsed = time.perf_counter() - time_init
    print(f'elapsed time: {time_elapsed:.3f} s | {time_elapsed / 12000:.5f} (s/call)')

    pd.DataFrame.from_records(d_calls).to_csv('eur_call_tau.csv')
