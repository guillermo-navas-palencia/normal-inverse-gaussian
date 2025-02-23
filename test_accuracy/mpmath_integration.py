from mpmath import *


def arg_mpmathify(x, alpha, beta, mu, delta):
    x = mpmathify(x)
    alpha = mpmathify(alpha)
    beta = mpmathify(beta)
    mu = mpmathify(mu)
    delta = mpmathify(delta)
    
    return x, alpha, beta, mu, delta


def normcdf(x):
    return erfc(-x / sqrt(mp.mpf('2'))) / mp.mpf('2')


def fun_phi(t, x, alpha, beta, mu, delta):
    two = mp.mpf('2')
    gamma = sqrt(alpha ** two - beta ** two)
    C = delta / sqrt(two * pi)
    return C * normcdf((x - (mu + beta * t)) / sqrt(t)) * t ** mp.mpf('-3/2') * exp(-(delta - gamma * t) ** two / two / t)


def quad_phi(x, alpha, beta, mu, delta, a=mp.zero, b=inf, digits=50):
    prec = mp.dps
    mp.dps = digits

    x, alpha, beta, mu, delta = arg_mpmathify(x, alpha, beta, mu, delta)

    cdf = quad(lambda t: fun_phi(t, x, alpha, beta, mu, delta), [a, b])

    mp.dps = prec

    return cdf