import numpy as np

from scipy import special


def _solve_expx_x_logx(tau, tol, max_steps=10):
    """Solves the equation

    log(pi/tau) = pi/2 * exp(x) - x - log(x)

    approximately using Newton's method. The approximate solution is guaranteed
    to overestimate.
    """
    # Initial guess
    x = np.log(2 / np.pi * np.log(np.pi / tau))
    # x = ln(tau / np.pi) - lambertw(-tau / 2, -1)
    # x = (
    #     mp.mpf(1) / 2
    #     - mp.log(mp.sqrt(mp.pi / tau))
    #     - mp.lambertw(-mp.sqrt(mp.exp(1) * mp.pi * tau) / 4, -1)
    # )

    def f0(x):
        return np.pi / 2 * np.exp(x) - x - np.log(x * np.pi / tau)

    def f1(x):
        return np.pi / 2 * np.exp(x) - 1 - 1 / x

    f0x = f0(x)
    success = False
    # At least one step is performed. This is required for the guarantee of
    # overestimation.
    for _ in range(max_steps):
        x -= f0x / f1(x)
        f0x = f0(x)
        if abs(f0x) < tol:
            success = True
            break

    assert success
    return x


def nig_integration(f, N, eps, max_steps=10):
    alpha2 = N / 2

    h = _solve_expx_x_logx(tau=eps ** 2, tol=1e-10)
    print(f'h={h}')

    last_error_estimate = None
    last_int = 0

    # Start iteration
    for level in range(max_steps + 1):
        print("lambert ", special.lambertw(-(eps ** 2) / h / 2, -1).real)
        j = int(np.log(-2 / np.pi * special.lambertw(-(eps ** 2) / h / 2, -1).real) / h)
        print(f'\nlevel: {level}. j: {j}')

        if level == 0:
            t = [0]
        else:
            t = h * np.arange(1, j + 1, 2)        

        # print(np.arange(1, j + 1, 2))
        sinh_t = np.pi / 2 * np.sinh(t)
        cosh_t = np.pi / 2 * np.cosh(t)
        cosh_sinh_t = np.cosh(sinh_t)
        exp_sinh_t = np.exp(sinh_t)

        y0 = alpha2 / exp_sinh_t / cosh_sinh_t
        y1 = -alpha2 * cosh_t / cosh_sinh_t ** 2
        weights = -h * y1

        # print(y0, y1, weights)

        fly = f(y0)
        fry = f(N - y0)

        # print(fly, fry)
        lsm = fly * weights
        rsm = fry * weights

        if level == 0:
            # The root level only contains one node, the midpoint; function values of
            # f_left and f_right are equal here. Deliberately take lsummands here.
            value_estimates = list(lsm)
        else:
            value_estimates.append(
                # Take the estimation from the previous step and half the step size.
                # Fill the gaps with the sum of the values of the current step.
                value_estimates[-1] / 2
                + np.sum(lsm)
                + np.sum(rsm)
            )        

        print(value_estimates)

        if abs(1 - last_int / value_estimates[-1]) < eps:
            success = True
            break

        last_int = value_estimates[-1]
        h /= 2

    return value_estimates[-1]
