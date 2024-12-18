import json
import os.path
import time

from ctypes import c_double
from typing import Dict, Tuple, Union

import numpy as np
import numpy.ctypeslib as npct
import numpy.typing as npt
import pandas as pd

from scipy import stats

from mpmath_nig import quad_nig


# load library
libabspath = os.path.dirname(os.path.abspath(__file__))

lib = npct.load_library('../../_nig.so', libabspath)

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


def run_test_set(
    name: str,
    x: npt.NDArray,
    alpha: npt.NDArray,
    beta: npt.NDArray,
    mu: npt.NDArray,
    delta: npt.NDArray
) -> None:
    
    filename_results = f'../results/{name}.csv'
    filename_data = f'../data/{name}.csv'

    with open(filename_results, 'w') as f:
        header = 'x,alpha,beta,mu,delta,mpmath,scipy,gnp\n'
        f.write(header)

    with open(filename_data, 'w') as f:
        header = 'x,alpha,beta,mu,delta,mpmath\n'
        f.write(header)

    n_samples = len(x)
    for i in range(n_samples):
        _x = x[i]
        _alpha = alpha[i]
        _beta = beta[i]
        _mu = mu[i]
        _delta = delta[i]

        # mpmath
        mp_result = quad_nig(
            x=_x,
            alpha=_alpha,
            beta=_beta,
            mu=_mu,
            delta=_delta,
            digits=100)

        # SciPy
        scipy_result = stats.norminvgauss.cdf(
            x=_x,
            a=_alpha*_delta,
            b=_beta*_delta,
            loc=_mu,
            scale=_delta)
        
        # C++ implementation
        gnp_result = nig_cdf(
            x=_x,
            alpha=_alpha,
            beta=_beta,
            mu=_mu,
            delta=_delta)
        
        instance_results = f'{_x},{_alpha},{_beta},{_mu},{_delta},{float(mp_result):.16E},{scipy_result:.16E},{gnp_result:.16E}\n'
        instance_data = f'{_x},{_alpha},{_beta},{_mu},{_delta},{float(mp_result):.16E}\n'

        with open(filename_results, 'a') as f:
            f.write(instance_results)

        with open(filename_data, 'a') as f:
            f.write(instance_data)

        if i % 100 == 0:
            print(f'{i} / {n_samples}')


def run_test_set_with_benchmark(name: str) -> None:
    filename_results = f'../results/{name}.csv'
    filename_data = f'../data/{name}.csv'

    df = pd.read_csv(filepath_or_buffer=filename_data)
    n_samples = df.shape[0]

    with open(filename_results, 'w') as f:
        header = 'x,alpha,beta,mu,delta,mpmath,scipy,gnp\n'
        f.write(header)

    cols = ['x', 'alpha', 'beta', 'mu', 'delta', 'mpmath']
    for i, row in df.iterrows():
        x, alpha, beta, mu, delta, mp_result = row[cols]

        scipy_result = stats.norminvgauss.cdf(
            x=x,
            a=alpha*delta,
            b=beta*delta,
            loc=mu,
            scale=delta)

        gnp_result = nig_cdf(
            x=x,
            alpha=alpha,
            beta=beta,
            mu=mu,
            delta=delta)

        instance_results = f'{x},{alpha},{beta},{mu},{delta},{mp_result},{scipy_result:.16E},{gnp_result:.16E}\n'

        with open(filename_results, 'a') as f:
            f.write(instance_results)

        if i % 100 == 0:
            print(f'{i} / {n_samples}')


def rerun_test_set_mpmath_accurate(name: str, max_err: float = 5e-13) -> None:
    filename_accurate = f'../results/{name}_mpmath_accurate.csv'
    filename_data = f'../data/{name}_mpmath_accurate.csv'

    # Read csv file
    df = pd.read_csv(filepath_or_buffer=f'../results/{name}.csv')
    
    # Replace NaN str with NaN
    df.loc[df['scipy'] == 'NAN', 'scipy'] = np.nan
    df['scipy'] = df['scipy'].astype(float)

    df.loc[df['gnp'] == 'NAN', 'gnp'] = np.nan
    df['gnp'] = df['gnp'].astype(float)    

    # Compute relative errors wrt mpmath
    df['relerr_gnp'] = np.absolute(df['gnp'] / df['mpmath'] - 1)
    df['abs_gnp'] = np.absolute(df['gnp'] - df['mpmath'])

    mask_inf = df['relerr_gnp'] == np.inf
    mask_nan = (df['relerr_gnp'].isna() & (df['mpmath'] != df['gnp']) & (df['mpmath'] != 0))
    mask_err = (df['relerr_gnp'] > max_err) | mask_nan | mask_inf

    with open(filename_accurate, 'w') as f:
        header = 'x,alpha,beta,mu,delta,mpmath,scipy,gnp\n'
        f.write(header)    
    
    with open(filename_data, 'w') as f:
        header = 'x,alpha,beta,mu,delta,mpmath\n'
        f.write(header)

    df_err = df[mask_err]
    n_samples = df_err.shape[0]

    cols = ['x', 'alpha', 'beta', 'mu', 'delta', 'mpmath']
    for it, (i, row) in enumerate(df_err.iterrows()):
        x, alpha, beta, mu, delta, mp_result = row[cols]

        # mpmath
        mp_result = quad_nig(
            x=x,
            alpha=alpha,
            beta=beta,
            mu=mu,
            delta=delta,
            digits=300)

        # SciPy
        scipy_result = stats.norminvgauss.cdf(
            x=x,
            a=alpha*delta,
            b=beta*delta,
            loc=mu,
            scale=delta)
        
        # C++ implementation
        gnp_result = nig_cdf(
            x=x,
            alpha=alpha,
            beta=beta,
            mu=mu,
            delta=delta)
        
        instance_results = f'{x},{alpha},{beta},{mu},{delta},{float(mp_result):.16E},{scipy_result:.16E},{gnp_result:.16E}\n'
        instance_data = f'{x},{alpha},{beta},{mu},{delta},{float(mp_result):.16E}\n'

        with open(filename_accurate, 'a') as f:
            f.write(instance_results)

        with open(filename_data, 'a') as f:
            f.write(instance_data)

        if it % 20 == 0:
            print(f'{it} / {n_samples}')


def create_test_set(
    n_samples: int,
    x_min: float,
    x_max: float,
    alpha_min: float,
    alpha_max: float,
    beta_min: float,
    beta_max: float,
    mu_min: float,
    mu_max: float,
    delta_min: float,
    delta_max: float,
    seed: int = 1
) -> Tuple[npt.NDArray]:
    rs = np.random.default_rng(seed=seed)

    if x_min == x_max:
        x = np.full(n_samples, x_min)
    else:
        x = rs.uniform(x_min, x_max, n_samples)

    alpha = rs.uniform(alpha_min, alpha_max, n_samples)

    if beta_min == beta_max:
        beta = np.full(n_samples, beta_min)
    else:
        ratio = rs.uniform(-0.999, 0.999, n_samples)
        beta = alpha * ratio
    
    if mu_min == mu_max:
        mu = np.full(n_samples, mu_min)
    else:
        mu = rs.uniform(mu_min, mu_max, n_samples)
    
    delta = rs.uniform(delta_min, delta_max, n_samples)

    return x, alpha, beta, mu, delta


def test_accuracy(name: str) -> None:
    with open(f'../data/{name}.json', 'r') as f:
        data = json.load(f)

    x, alpha, beta, mu, delta = create_test_set(
        n_samples=data['n_samples'],
        x_min=data['x_min'],
        x_max=data['x_max'],
        alpha_min=data['alpha_min'],
        alpha_max=data['alpha_max'],
        beta_min=data['beta_min'],
        beta_max=data['beta_max'],
        mu_min=data['mu_min'],
        mu_max=data['mu_max'],
        delta_min=data['delta_min'],
        delta_max=data['delta_max']
    )

    run_test_set(
        name=name,
        x=x,
        alpha=alpha,
        beta=beta,
        mu=mu,
        delta=delta
    )


def generate_accuracy_summary_beta_eq_zero(
    name: str,
    max_err: float = 5e-13
) -> None:
    # Read csv file
    df = pd.read_csv(filepath_or_buffer=f'../results/{name}.csv')
    n_samples = df.shape[0]

    # Replace NaN str with NaN
    df.loc[df['scipy'] == 'NAN', 'scipy'] = np.nan
    df['scipy'] = df['scipy'].astype(float)

    df.loc[df['gnp'] == 'NAN', 'gnp'] = np.nan
    df['gnp'] = df['gnp'].astype(float)

    # Compute relative errors wrt mpmath
    df['relerr_gnp'] = np.absolute(df['gnp'] / df['mpmath'] - 1)
    df['abs_gnp'] = np.absolute(df['gnp'] - df['mpmath'])

    df['relerr_scipy'] = np.absolute(df['scipy'] / df['mpmath'] - 1)
    df['abs_scipy'] = np.absolute(df['scipy'] - df['mpmath'])

    # scipy error summary
    mask_inf = df['relerr_scipy'] == np.inf
    mask_nan = (df['relerr_scipy'].isna() & (df['mpmath'] != df['scipy']) & (df['mpmath'] != 0))
    mask_err = (df['relerr_scipy'] > max_err) | mask_nan
    
    print('\nSUMMARY SciPy')
    print('-------------')
    print(f'Cases inf: {mask_inf.sum()}')
    print(f'Cases NaN: {mask_nan.sum()}')

    n_errors = np.count_nonzero(~mask_err)
    print(f'Cases err < {max_err}: {n_errors} ({n_errors / n_samples:.2%})')

    # gnp error summary
    mask_inf = df['relerr_gnp'] == np.inf
    mask_nan = (df['relerr_gnp'].isna() & (df['mpmath'] != df['gnp']) & (df['mpmath'] != 0))
    mask_err = (df['relerr_gnp'] > max_err) | mask_nan

    print('SUMMARY gnp')
    print('-----------')
    print(f'Cases inf: {mask_inf.sum()}')
    print(f'Cases NaN: {mask_nan.sum()}')
    
    n_errors = np.count_nonzero(~mask_err)
    print(f'Cases err < {max_err}: {n_errors} ({n_errors / n_samples:.2%})')

    # Algorithmic choice
    df['target'] = (~mask_err).astype(int)

    df['omega'] = np.sqrt(df['delta'] ** 2 + (df['x'] - df['mu']) ** 2)
    df['aw'] = df['alpha'] * df['omega']
    df['xmu'] = df['x'] - df['mu']
    df['xmu2'] = df['xmu'] ** 2
    df['aow'] = df['alpha'] / df['omega']
    df['woa'] = df['omega'] / df['alpha']
    df['da'] = df['delta'] * df['alpha']

    df['method'] = 'integration'

    criteria_asymp_xmu = (df['xmu2'] >= 70) & (df['aow'] >= 1.0)
    df.loc[criteria_asymp_xmu, 'method'] = 'asymptotic_xmu'

    criteria_asymp_alpha_c1 = (df['xmu2'] <= 2.5) & (df['da'] >= 200)
    criteria_asymp_alpha_c2 = (df['alpha'] >= 5) & (df['delta'] >= 10)
    criteria_asymp_alpha = criteria_asymp_alpha_c1 & criteria_asymp_alpha_c2
    df.loc[criteria_asymp_alpha, 'method'] = 'asymptotic_alpha'

    criteria1 = (df['xmu'].abs() <= 20) & (df['aow'] <= 0.25) & (df['delta'] >= 2.0 * df['xmu'].abs())
    criteria2 = (df['xmu2'] <= 1.25) & (df['aw'] <= 750)
    criteria3 = (df['delta'] >= 1.0)

    df.loc[(criteria1 | criteria2) & criteria3, 'method'] = 'series'

    print(df['method'].value_counts())
    print(df.groupby('method')[['target', 'relerr_gnp']].describe().to_string())


def generate_accuracy_summary_x_eq_mu(
    name: str,
    max_err: float = 5e-13
) -> None:
    # Read csv file
    df = pd.read_csv(filepath_or_buffer=f'../results/{name}.csv')
    n_samples = df.shape[0]

    # Replace NaN str with NaN
    df.loc[df['scipy'] == 'NAN', 'scipy'] = np.nan
    df['scipy'] = df['scipy'].astype(float)

    df.loc[df['gnp'] == 'NAN', 'gnp'] = np.nan
    df['gnp'] = df['gnp'].astype(float)

    # Compute relative errors wrt mpmath
    df['relerr_gnp'] = np.absolute(df['gnp'] / df['mpmath'] - 1)
    df['abs_gnp'] = np.absolute(df['gnp'] - df['mpmath'])

    df['relerr_scipy'] = np.absolute(df['scipy'] / df['mpmath'] - 1)
    df['abs_scipy'] = np.absolute(df['scipy'] - df['mpmath'])

    # scipy error summary
    mask_inf = df['relerr_scipy'] == np.inf
    mask_nan = (df['relerr_scipy'].isna() & (df['mpmath'] != df['scipy']) & (df['mpmath'] != 0))
    mask_err = (df['relerr_scipy'] > max_err) | mask_nan
    
    print('\nSUMMARY SciPy')
    print('-------------')
    print(f'Cases inf: {mask_inf.sum()}')
    print(f'Cases NaN: {mask_nan.sum()}')

    n_errors = np.count_nonzero(~mask_err)
    print(f'Cases err < {max_err}: {n_errors} ({n_errors / n_samples:.2%})')

    # gnp error summary
    mask_inf = df['relerr_gnp'] == np.inf
    mask_nan = (df['relerr_gnp'].isna() & (df['mpmath'] != df['gnp']) & (df['mpmath'] != 0))
    mask_err = (df['relerr_gnp'] > max_err) | mask_nan

    print('SUMMARY gnp')
    print('-----------')
    print(f'Cases inf: {mask_inf.sum()}')
    print(f'Cases NaN: {mask_nan.sum()}')
    
    n_errors = np.count_nonzero(~mask_err)
    print(f'Cases err < {max_err}: {n_errors} ({n_errors / n_samples:.2%})')

    # Algorithmic choice
    df['target'] = (~mask_err).astype(int)
    df['rba'] = df['beta'] / df['alpha']
    df['method'] = 'integration'

    criteria1 = (df['alpha'] <= 10.0) & (df['delta'] <= 10.0)
    criteria2 = (df['beta'] <= 1.5) & (df['rba'] <= 0.9)
    df.loc[criteria1 & criteria2, 'method'] = 'series'

    print(df['method'].value_counts())
    print(df.groupby('method')[['target', 'relerr_gnp']].describe().to_string())


def run_test_set_timing(name: str) -> None:
    filename_results = f'../results/{name}.csv'
    df = pd.read_csv(filepath_or_buffer=filename_results)

    x = df['x'].values
    alpha = df['alpha'].values
    beta = df['beta'].values
    mu = df['mu'].values
    delta = df['delta'].values

    time_init = time.perf_counter()
    for i in range(df.shape[0]):
        stats.norminvgauss.cdf(
            x=x[i],
            a=alpha[i]*delta[i],
            b=beta[i]*delta[i],
            loc=mu[i],
            scale=delta[i])
    elapsed_scipy = time.perf_counter() - time_init

    time_init = time.perf_counter()
    for i in range(df.shape[0]):
        nig_cdf(
            x=x[i],
            alpha=alpha[i],
            beta=beta[i],
            mu=mu[i],
            delta=delta[i])

    elapsed_cpp = time.perf_counter() - time_init

    print(f'scipy: {elapsed_scipy:.4f}s')
    print(f'cpp  : {elapsed_cpp:.4f}s')

if __name__ == '__main__':
    name = 'test_x_eq_mu_large'

    # 1. Test accuracy comparing with mpmath and SciPy implementations
    # test_accuracy(name=name)
    run_test_set_with_benchmark(name=name)
    # rerun_test_set_mpmath_accurate(name=name)

    # 2. Generate test summary
    # generate_accuracy_summary_beta_eq_zero(name=name)
    # generate_accuracy_summary_x_eq_mu(name=name)

    # 3. Timing SciPy vs C++ via ctypes
    # run_test_set_timing(name=name)

