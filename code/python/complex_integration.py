import numpy as np

from scipy import integrate, special


def quad_complex(x, alpha, beta, mu, delta, n):
    gamma = np.sqrt(alpha ** 2 - beta ** 2)
    da = delta * alpha

    C = 1 / np.pi

    xk, wk = special.roots_laguerre(n=n)
    print(xk, wk)

    qk = da + xk
    rk = qk ** 2 / delta ** 2 - alpha ** 2
    vk = np.sqrt(rk) * 1j
    sk = vk - beta

    integral = wk @ (np.exp(-(x-mu) * sk + delta * gamma - da -np.log(sk) - np.log(vk) + np.log(qk)))

    return 1.0 + C / delta ** 2 * np.imag(integral)



def quad_complex_direct(x, alpha, beta, mu, delta):
    gamma = np.sqrt(alpha ** 2 - beta ** 2)
    da = delta * alpha
    xmu = x - mu

    C = 1 / np.pi

    def f(q):
        r = q ** 2 / delta ** 2 - alpha ** 2
        v = np.sqrt(r) * 1j
        s = v - beta
        return np.exp(-xmu * s - q + delta * gamma -np.log(s) - np.log(v) + np.log(q))
    
    integral = integrate.quad(func=f, a=da, b=np.inf, complex_func=True, full_output=1, epsabs=1e-14)
    # print(integral[1])

    return 1.0 + C / delta ** 2 * np.imag(integral[0])