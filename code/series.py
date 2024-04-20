import numpy as np

from scipy import special


def series_x_mu_b_zero(
    x: float,
    alpha: float,
    mu: float,
    delta: float,
    maxiter: int = 200,
    eps: float = 5e-15
) -> tuple[float, int]:
    
    # Parameters
    xmu = x - mu
    xmu2 = xmu * xmu
    ad = alpha * delta
    aod = alpha / delta

    # Constant
    C = delta * np.exp(delta * alpha) / np.pi

    # Series: compute the first two terms with K_1 and K_2 for recursion
    kn = special.kn(1, ad)
    knp1 = special.kn(2, ad)

    num = xmu * aod
    s = num * kn

    num *= -xmu2 * aod / 2.0
    s += num * knp1 / 3.0

    sp = s
    for k in range(2, maxiter):
        # Bessel recursion
        knn = kn + 2 * k / ad * knp1

        # New term
        num *= -xmu2 * aod / (2 * k)
        s += num * knn / (2 * k + 1)

        # Check convergence
        if abs(1.0 - sp / s) < eps:
            return 0.5 + C * s, k
        else:
            kn = knp1
            knp1 = knn
            sp = s

    # No convergence
    return -1, maxiter


def series_x_mu_b_zero_pos(
    x: float,
    alpha: float,
    mu: float,
    delta: float,
    maxiter: int = 200,
    eps: float = 5e-15
) -> tuple[float, int]:

    # Parameters
    xmu = x - mu
    xmu2 = xmu * xmu

    # Use hypotenuse
    omega = np.sqrt(xmu2 + delta*delta)
    aw = alpha * omega
    aow = alpha / omega

    # Constant
    C = delta * np.exp(delta * alpha) / np.pi

    # Series: compute the first two terms with K_1 and K_2 for recursion
    kn = special.kn(1, aw)
    knp1 = special.kn(2, aw)

    num = xmu * aow
    s = num * kn

    num *= xmu2 * aow / 3
    s += num * knp1

    sp = s
    for k in range(2, maxiter):
        # Bessel recursion
        knn = kn + 2 * k / aw * knp1

        # New term
        num *= xmu2 * aow / (2 * k + 1)
        s += num * knn

        # Check convergence
        if abs(1.0 - sp / s) < eps:
            return 0.5 + C * s, k
        else:
            kn = knp1
            knp1 = knn
            sp = s        

    # No convergence
    return -1, maxiter


def asymtotic_x_mu_b_zero(
    x: float,
    alpha: float,
    mu: float,
    delta: float,
    maxiter: int = 200,
    eps: float = 5e-15
) -> tuple[float, int]:

    # Parameters
    xmu = x - mu
    xmu2 = xmu * xmu

    # Use hypotenuse
    omega = np.sqrt(xmu2 + delta*delta)
    aw = alpha * omega
    woa = omega / alpha

    # Constant
    C = delta * np.exp(delta * alpha) / np.pi

    # Series: compute the first two terms with K_0 and K_1
    kn = special.kn(0, aw)
    knp1 = special.kn(1, aw)

    oxmu2 = 1.0 / xmu2
    num = -1.0 / xmu
    s = num * kn

    num *= -woa * oxmu2
    s += num * knp1

    sp = s
    for k in range(2, maxiter):
        # Bessel recursion
        knn = kn + 2 * (k - 1) / aw * knp1

        # New term
        num *= -woa * oxmu2 * (2 * k - 1)
        s += num * knn

        # Check convergence
        if abs(1.0 - sp / s) < eps:
            return C * s, k
        else:
            kn = knp1
            knp1 = knn
            sp = s

    # No convergence
    return -1, maxiter        