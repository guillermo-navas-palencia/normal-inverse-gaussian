import json
import os.path
import time

from ctypes import c_double

import numpy as np
import numpy.ctypeslib as npct
import numpy.typing as npt

from scipy import special, stats

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


def analytical_digital_option_cash_or_nothing(
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

    cdf = nig_cdf(x=k, alpha=alpha, beta=-beta, mu=0.0, delta=delta * tau)

    call = np.exp(-r * tau) * cdf
    put = np.exp(-r * tau) * (1.0 - cdf)

    return call, put


def analytical_digital_option_asset_or_nothing(
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

    cdf = nig_cdf(x=k, alpha=alpha, beta=-beta-1, mu=0.0, delta=delta * tau)
    call = St * np.exp(-q * tau) * cdf
    put = St * np.exp(-q * tau) - call

    return call, put


def analytical_european_option(
    St: float,
    K: float,
    r: float,
    q: float,
    tau: float,
    alpha: float,
    beta: float,
    delta: float
) -> float:
    digital_cn_call, _ = analytical_digital_option_cash_or_nothing(
        St=St, K=K, r=r, q=q, tau=tau, alpha=alpha, beta=beta, delta=delta
    )

    # digital option asset-or-nothing (tau=1)
    digital_an_call, _ = analytical_digital_option_asset_or_nothing(
        St=St, K=K, r=r, q=q, tau=tau, alpha=alpha, beta=beta, delta=delta
    )

    call = digital_an_call - K * digital_cn_call
    put = call - St + K * np.exp(-r * tau)

    return call, put


def analytical_log_option(
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

    cdf_nig = nig_cdf(x=k, alpha=alpha, beta=-beta, mu=0.0, delta=delta * tau)

    gh_dist = stats.genhyperbolic(
        p=0.5,
        a=alpha * delta * tau,
        b=-beta * delta * tau,
        loc=0.0,
        scale=delta * tau)

    cdf_gh = gh_dist.cdf(k)

    H1 = (delta * tau * np.exp(gamma * delta * tau - beta * k)) / np.pi
    H2 = beta * delta * tau / gamma

    call = np.exp(-r * tau) * (
        k * cdf_nig +
        H1 * special.k0(alpha * np.sqrt(k**2 + (delta*tau)**2)) +
        H2 * cdf_gh
    )

    put = call - np.exp(-r * tau) * k

    return call, put


def analytical_power_option_calls(
    St: float,
    K: float,
    r: float,
    q: float,
    a: float,
    tau: float,
    alpha: float,
    beta: float,
    delta: float
) -> tuple[float, ...]:
    gamma = np.sqrt(alpha**2 - beta**2)
    omega = delta * (np.sqrt(alpha**2 - (beta + 1)**2) - gamma)
    ka = np.log(St / K ** (1/a)) + (r - q + omega) * tau
    discount = np.exp(-r * tau)

    cdf = nig_cdf(x=ka, alpha=alpha, beta=-beta, mu=0.0, delta=delta * tau)
    cn_call = discount * cdf
    cn_put = discount * (1.0 - cdf)

    wa = delta * tau * (np.sqrt(alpha**2 - (beta + a)**2) - gamma)
    an_call = K * np.exp(-r * tau + a * ka - wa) * nig_cdf(x=ka, alpha=alpha, beta=-beta-a, mu=0.0, delta=delta * tau)
    an_put = 0.0

    european_call = an_call - K * cn_call
    european_put = 0.0

    return cn_call, cn_put, an_call, an_put, european_call, european_put


def mc(
    St: float,
    K: float,
    r: float,
    q: float,
    a: float,
    tau: float,
    alpha: float,
    beta: float,
    delta: float,
    n_paths: int = 1_000_000,
    seed: int | None = None
) -> float:
    
    nig_dist = stats.norminvgauss(
        a=alpha * delta * tau,
        b=beta * delta * tau,
        loc=0.0,
        scale=delta * tau
    )

    Xt = nig_dist.rvs(size=n_paths, random_state=seed)
    gamma = np.sqrt(alpha**2 - beta**2)
    omega = delta * (np.sqrt(alpha**2 - (beta + 1)**2) - gamma)

    ST = St * np.exp((r - q + omega) * tau + Xt)

    # digital option cash-or-nothing
    discount = np.exp(-r * tau)
    digital_cn_call = discount * (ST > K).mean()
    digital_cn_put = discount * (ST < K).mean()

    # digital option asset-or-nothing
    digital_an_call = discount * (ST * (ST > K)).mean()
    digital_an_put = discount * (ST * (ST < K)).mean()

    # european option
    european_call = discount * np.maximum(ST - K, 0).mean()
    european_put = discount * np.maximum(K - ST, 0).mean()

    # log option
    log_call = discount * np.maximum(np.log(ST) - np.log(K), 0).mean()
    log_put = discount * np.maximum(np.log(K) - np.log(ST), 0).mean()

    # power option
    power_cn_call = discount * (ST ** a > K).mean()
    power_cn_put = discount * (ST ** a < K).mean()
    power_an_call = discount * (ST ** a * (ST ** a > K)).mean()
    power_an_put = discount * (ST ** a * (ST ** a < K)).mean()
    power_european_call = discount * np.maximum(ST ** a - K, 0).mean()
    power_european_put = discount * np.maximum(K - ST ** a, 0).mean()

    print('\ndigital option')
    print(f'cash-or-nothing (CALL): {digital_cn_call:.4f}')
    print(f'cash-or-nothing (PUT): {digital_cn_put:.4f}')
    print(f'asset-or-nothing (CALL): {digital_an_call:.4f}')
    print(f'asset-or-nothing (PUT): {digital_an_put:.4f}')

    print('\neuropean option')
    print(f'CALL: {european_call:.4f}')
    print(f'PUT: {european_put:.4f}')

    print('\nlog option')
    print(f'CALL: {log_call:.4f}')
    print(f'PUT: {log_put:.4f}')    

    print('\npower option')
    print(f'cash-or-nothing (CALL): {power_cn_call:.4f}')
    print(f'cash-or-nothing (PUT): {power_cn_put:.4f}')
    print(f'asset-or-nothing (CALL): {power_an_call:.4f}')
    print(f'asset-or-nothing (PUT): {power_an_put:.4f}')
    print(f'european (CALL): {power_european_call:.4f}')
    print(f'european (PUT): {power_european_put:.4f}')


if __name__ == '__main__':
    St = 3000
    K = 4000
    r = 0.01
    q = 0.0
    a = 1.2

    alpha = 8.9932
    beta = -4.5176
    delta = 1.1528

    # Monte Carlo
    print('\nMonte Carlo')
    print('-------------------')
    mc(St=St, K=K, r=r, q=q, a=a, tau=2, alpha=alpha, beta=beta, delta=delta,
       n_paths=10_000_000, seed=1)
    
    # Analytical
    print('\nAnalytical')
    print('-------------------')

    # digital option cash-or-nothing (tau=2)
    digital_cn_call, digital_cn_put = analytical_digital_option_cash_or_nothing(
        St=St, K=K, r=r, q=q, tau=2, alpha=alpha, beta=beta, delta=delta
    )

    # digital option asset-or-nothing (tau=1)
    digital_an_call, digital_an_put = analytical_digital_option_asset_or_nothing(
        St=St, K=K, r=r, q=q, tau=1, alpha=alpha, beta=beta, delta=delta
    )

    # vanilla european option
    european_call, european_put = analytical_european_option(
        St=St, K=K, r=r, q=q, tau=1, alpha=alpha, beta=beta, delta=delta
    )

    # log option
    log_call, log_put = analytical_log_option(
        St=St, K=K, r=r, q=q, tau=2, alpha=alpha, beta=beta, delta=delta
    )

    # power option
    (power_cn_call, power_cn_put, power_an_call, power_an_put,
     power_european_call, power_european_put) = analytical_power_option_calls(
        St=St, K=K, r=r, q=q, a=a, tau=2, alpha=alpha, beta=beta, delta=delta
    )

    print('\ndigital')
    print(f'cash-or-nothing (CALL): {digital_cn_call:.8f}')
    print(f'cash-or-nothing (PUT): {digital_cn_put:.8f}')
    print(f'asset-or-nothing (CALL): {digital_an_call:.8f}')
    print(f'asset-or-nothing (PUT): {digital_an_put:.8f}')

    print('\neuropean option')
    print(f'CALL: {european_call:.8f}')
    print(f'PUT: {european_put:.8f}')

    print('\nlog option')
    print(f'CALL: {log_call:.8f}')
    print(f'PUT: {log_put:.8f}')

    print('\npower option')
    print(f'cash-or-nothing (CALL): {power_cn_call:.8f}')
    print(f'cash-or-nothing (PUT): {power_cn_put:.8f}')
    print(f'asset-or-nothing (CALL): {power_an_call:.8f}')
    print(f'asset-or-nothing (PUT): {power_an_put:.8f}')
    print(f'european (CALL): {power_european_call:.8f}')
    print(f'european (PUT): {power_european_put:.8f}')