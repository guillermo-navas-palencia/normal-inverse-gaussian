import numpy as np

from scipy import special


def series_x_mu_b0(
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

    num = C * xmu * aod
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
            return 0.5 + s, k
        else:
            kn = knp1
            knp1 = knn
            sp = s

    # No convergence
    return -1, maxiter


def series_x_mu_b0_pos(
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

    num = C * xmu * aow
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
        # print(k, s, num, knn)

        # Check convergence
        if abs(1.0 - sp / s) < eps:
            return 0.5 + s, k
        else:
            kn = knp1
            knp1 = knn
            sp = s        

    # No convergence
    return -1, maxiter


def series_x_mu_b0_pos_robust(
    x: float,
    alpha: float,
    mu: float,
    delta: float,
    maxiter: int = 1000,
    eps: float = 5e-15    
) -> tuple[float, int]:

    # Parameters
    xmu = x - mu
    xmu2 = xmu * xmu

    # Use hypotenuse
    omega = np.sqrt(xmu2 + delta*delta)
    aw = alpha * omega
    aow = alpha / omega

    z = xmu2 * aow
    z2 = z*z
    awz = aw * z
    awz2 = awz * z
    awt5 = aw * 5
    awt2 = aw * 2
    zt4 = z * 4
    zt2 = z * 2

    # Constant
    C = delta * np.exp(delta * alpha) / np.pi * xmu * aow

    # Series: compute first two partial sums
    y0 = 0
    y1 = special.kn(1, aw)
    y2 = y1 + special.kn(2, aw) * z / 3

    sp = y2
    for k in range(maxiter):
        # Improve numerical stability
        # n = (3 + 2*k)
        # num = -awz2 * y0 + z*(-12 -14*k -4*k*k + awz) * y1 + n * (awt5 + k*awt2 + zt4 + zt2*k) * y2
        # den = n * (5 + 2 * k) * aw

        # num = -awz2 * y0 + z* (-2 * (k+2)* (3 + 2*k) + awz) * y1 + (3 + 2*k) * ((5 + 2*k) * aw + 2*(2+k)*z) * y2
        # den = (3 + 2*k) * (5 + 2*k) * aw
        # s = num / den

        n0 = -z2 / ((3 + 2*k) * (5 + 2*k)) * y0
        # n1 = z* (-2 * (k+2)* (3 + 2*k) + awz) * y1 / ((3 + 2*k) * (5 + 2*k) * aw)

        n11 = z* (-2 * (k+2)) / ((5 + 2*k) * aw)
        n12 = z2 / ((3 + 2*k) * (5 + 2*k))
        n1 = (n11 + n12) * y1

        n2 = (1 + (2*(2+k)*z) / ((5 + 2*k) * aw)) * y2 

        s = n0 + n1 + n2

        if abs(1.0 - sp / s) < eps:
            return 0.5 + C * s, k
        else:
            y0 = y1
            y1 = y2
            y2 = s
            sp = s

    # No convergence
    return -1, maxiter


def asymtotic_x_mu_b0_neg(
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
        # print(k, s)

        # Check convergence
        if abs(1.0 - sp / s) < eps:
            return C * s, k
        else:
            kn = knp1
            knp1 = knn
            sp = s

    # No convergence
    return -1, maxiter


def asymptotic_alpha_b0(
    x: float,
    alpha: float,
    mu: float,
    delta: float,
    maxiter: int = 200,
    eps: float = 5e-15
) -> tuple[float, int]:

    xmu = x - mu

    mirror = xmu < 0

    a = xmu / np.sqrt(2)
    b = delta / alpha

    z = alpha * alpha / 2
    h = -1 / 2
    xi = b * z

    # Constant
    C = delta * np.exp(delta * alpha) / np.sqrt(2 * np.pi) * z ** (-h)

    # Compute the first two terms of the Phi((x-mu) / sqrt(t)) at t=b
    c0 = special.erfc(-a / np.sqrt(b)) / 2
    c1 = -a/2 * np.exp(-a*a / b) / (np.sqrt(np.pi) * b ** (3/2))

    # Compute first elements of binomial sum of Bessel functions recursion
    kh = special.kv(h, 2 * xi)
    khp2 = special.kv(h + 2, 2 * xi)

    # print(kh)
    # print(khp2)

    q0 = 2 * xi ** h * kh
    q1 = 0
    q2 = 2 * xi ** (h + 2) * (khp2 - kh)

    s = c0 * q0
    ovz = 1 / z
    num = ovz

    print(1, s)
    # print(0, c0)
    # print(1, c1)

    sp = s
    for k in range(2, maxiter):
        # Phi recursion
        n = k - 2
        ck = ((n + 1) * c1 * (2*a**2 - 4*b*n - 3*b) - (2*n**2 + n) * c0) / (2*b**2 * (n+1) * (n+2))

        # ck2 = ((k - 1) * c1 * (xmu**2 - 4*b*(k - 2) - 3*b) - (2*(k - 2)**2 + k - 2) * c0) / (2*b**2 * (k - 1) * (k))
        # ck2 = (xmu**2 - 4*b*(k-2) - 3*b) / (2 * b**2 * k) * c1 - (2*(k - 2)**2 + k - 2) * c0 / (2*b**2 * (k - 1) * (k))
        # print(abs(ck))

        # Bessel recursion
        if k >= 3:
            n = k - 2
            qk = (n + h + 1 - 2 * xi) * q2 + xi * (2 * n + h + 1) * q1 + n * xi ** 2 * q0
        else:
            qk = q2

        # New term
        num *= ovz
        s += num * ck * qk
        # print(k, ck * qk * num)
        print(k, abs(num * ck * qk))

        if abs(1.0 - sp / s) < eps or k == maxiter - 1:
            return C * s, k
        else:
            c0 = c1
            c1 = ck

            if k >= 3:
                q0 = q1
                q1 = q2
                q2 = qk

            sp = s

    return -1, maxiter