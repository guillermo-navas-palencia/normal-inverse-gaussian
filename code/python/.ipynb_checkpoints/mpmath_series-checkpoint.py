from mpmath import *


def arg_mpmathify(x, alpha, mu, delta):
    x = mpmathify(x)
    alpha = mpmathify(alpha)
    mu = mpmathify(mu)
    delta = mpmathify(delta)
    
    return x, alpha, mu, delta


def normcdf(x):
    return erfc(-x / sqrt(2)) / 2


def fun_phi(t, x, alpha, mu, delta):
    C = delta / sqrt(2 * pi)
    return C * normcdf((x - mu) / sqrt(t)) * t ** (-3/2) * exp(-(delta - alpha * t) ** 2 / 2 / t)


def quad_phi(x, alpha, mu, delta, digits=50):
    prec = mp.dps
    mp.dps = digits
    
    x, alpha, mu, delta = arg_mpmathify(x, alpha, mu, delta)

    cdf = quad(lambda t: fun_phi(t, x, alpha, mu, delta), [0, inf])

    mp.dps = prec

    return cdf


def estimate_N_series_x_mu_b0_pos(x, alpha, mu, delta, epsilon=1e-16):
    omega = sqrt(delta**2 + (x-mu)**2)
    A = (x-mu) / omega
    B = A * delta * exp(delta * alpha) / (2 * pi * omega)
    C = epsilon**2 / (B**2 * pi)
    D = -4*log(A) / (C * A**2)

    L = log(D)
    approx = (-L + log(L)) /4/log(A) + 1/2

    print(float(A), float(B), float(C), float(D), float(approx))

    return mp.re(lambertw(-1, D)) / 4 / log(A) + 1/2


def estimate_N_series_x_mu_b0_pos_fast(x, alpha, mu, delta, epsilon=1e-16):
    omega = sqrt(delta**2 + (x-mu)**2)

    A = exp(1/2 - alpha * omega) * sqrt(pi / (2 * alpha * omega)) / sqrt(2)
    A *= (x - mu) * alpha / omega * delta * exp(delta * alpha) / pi
    A = epsilon / A

    B = (x - mu)**2 * alpha / omega

    return -log(A) / lambertw(-2 * log(A) / B / sqrt(e))


def estimate_Tk_series_x_mu_b0_pos_fast(x, alpha, mu, delta, N):
    omega = sqrt(delta**2 + (x-mu)**2)

    A = exp(- alpha * omega) * sqrt(pi / (2 * alpha * omega))
    A *= (x - mu) * alpha / omega * delta * exp(delta * alpha) / pi
    B = (x - mu)**2 * alpha / omega

    return A * B**N * exp(N + 1/2) / sqrt(2) / (2*N)**(N)
    # return -log(A) / lambertw(-2 * log(A) / B / sqrt(e))


def bound_series_x_mu_b0_pos(x, alpha, mu, delta, N):
    x, alpha, mu, delta = arg_mpmathify(x, alpha, mu, delta)
    omega = sqrt(delta**2 + (x-mu)**2)

    TN = 1/(2 * alpha * omega) * ((x - mu) / omega) ** (2 * N)
    TN *= sqrt(pi / (N + 1/2))

    C = ((x-mu) / omega) ** 2 / (2 * N + 3)
    C *= (N + 3/2 + sqrt((N + 3/2)**2 + (alpha * omega)**2))

    bound = TN / (1 - C)

    # Multiply normalizing factor
    bound *= (x - mu) * alpha / omega * delta * exp(delta * alpha) / pi

    return float(bound)


def bound_series_x_mu_b0_pos_alt(x, alpha, mu, delta, N):
    x, alpha, mu, delta = arg_mpmathify(x, alpha, mu, delta)
    omega = sqrt(delta**2 + (x-mu)**2)

    TN = 1/(2 * alpha * omega) * ((x - mu) / omega) ** (2 * N)
    TN *= sqrt(pi / (N + 1/2))

    C = 1 / (1 - ((x-mu) / omega)**2)
    bound = TN * C

    # Multiply normalizing factor
    bound *= (x - mu) * alpha / omega * delta * exp(delta * alpha) / pi

    return float(bound)


def series_x_mu_b0_pos(x, alpha, mu, delta, N, digits=50):
    prec = mp.dps
    mp.dps = digits

    x, alpha, mu, delta = arg_mpmathify(x, alpha, mu, delta)

    omega = sqrt(delta**2 + (x-mu)**2)
    z = (x-mu)**2 * alpha / omega

    s = mp.zero
    for k in range(N):
        r = z**k / fac2(2*k + 1)
        q = besselk(k + 1, alpha * omega)
        s += r * q

        # estimate1 = sqrt(pi / (2 * alpha * omega)) * (x-mu)**k * alpha ** 2 * exp(-alpha * omega) / fac2(2*k + 1)
        # estimate2 = sqrt(pi) / (2 * alpha * omega) * gamma(k + 1) / gamma(k + 3/2)
        # print(k, float(r * q), float(estimate1), float(estimate2))
        # estimate = sqrt(pi / e) * (alpha * omega)**(-k - 1) / sqrt(4*k + 2) * ((2*k + 2) / (2*k + 1)) ** (k + 1/2)
        # estimate1 = sqrt(pi / e) / sqrt(4*k + 2) * ((2*k + 2) / (2*k + 1)) ** (k + 1/2) * 1/(alpha * omega) * ((x-mu) / omega)**(2*k)
        # estimate2 = estimate * z**k
        # print(k, float(r*q), float(estimate1), float(estimate2), float(((x-mu) / omega)**(2*k)))

    s *= (x - mu) * alpha / omega
    C = delta * exp(delta * alpha) / pi
    cdf = s * C + mpf('1/2')

    mp.dps = prec

    return cdf


def estimate_N_series_x_mu_b0_alternating(x, alpha, mu, delta, epsilon=1e-16):
    A = (x - mu) / delta
    B = A * exp(delta * alpha) / pi
    C = epsilon / B

    D = -log(A) / (C * A)
    print(float(D))

    return (mp.re(lambertw(-1, D)) + log(A)) / (2 * log(A))


def bound_series_x_mu_b0_alternating(x, alpha, mu, delta, N):
    x, alpha, mu, delta = arg_mpmathify(x, alpha, mu, delta)

    TN = 1/(alpha * delta) * ((x - mu) / delta) ** (2 * N)
    TN *= 1 / (2*N + 1) 

    Y = 2 * (N + 1) * (2*N + 3) / (2*N + 1)

    C = ((x-mu) / delta) ** 2 / Y
    C *= (N + 3/2 + sqrt((N + 3/2)**2 + (alpha * delta)**2))

    bound = TN / (1 - C)

    # Multiply normalizing factor
    bound *= (x - mu) * alpha * exp(delta * alpha) / pi

    return float(bound)


def series_x_mu_b0_alternating(x, alpha, mu, delta, N, digits=50):
    prec = mp.dps
    mp.dps = digits

    x, alpha, mu, delta = arg_mpmathify(x, alpha, mu, delta)

    z = (x-mu)**2 * alpha / delta

    s = mp.zero
    for k in range(N):
        r = (-z)**k / (2**k * factorial(k) * (2*k + 1))
        q = besselk(k + 1, alpha * delta)
        s += r * q

    s *= (x - mu) * alpha / delta
    C = delta * exp(delta * alpha) / pi
    cdf = s * C + mpf('1/2')

    mp.dps = prec

    return cdf    