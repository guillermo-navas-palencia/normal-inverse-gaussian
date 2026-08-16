"""SciPy-compatible Normal Inverse Gaussian distribution.

Follows the same parametrisation as ``scipy.stats.norminvgauss`` so that
this package can be used as a drop-in replacement::

    # scipy
    from scipy import stats
    stats.norminvgauss.cdf(x, a=4.0, b=-0.8, loc=1.75, scale=2.0)

    # nig (identical call signature)
    from nig import norminvgauss
    norminvgauss.cdf(x, a=4.0, b=-0.8, loc=1.75, scale=2.0)

Parametrisation
---------------
SciPy uses shape parameters ``a`` and ``b`` together with the standard
``loc`` / ``scale`` shift-scale convention:

    a = alpha * delta   (tail heaviness,  a > |b|)
    b = beta  * delta   (asymmetry)
    loc   = mu          (location)
    scale = delta       (scale, > 0)

The Wikipedia / direct parametrisation (alpha, beta, mu, delta) is
available via :func:`nig.nig_cdf`.
"""

import numpy as np

from ._nig import nig_cdf as _nig_cdf_scalar


# Vectorised scalar kernel – handles Python scalars, lists, and NumPy arrays.
# otypes=[float] skips numpy's type-probing call (which would call the function
# with dummy inputs and could itself trigger the FPE warning).
_nig_cdf_vec = np.vectorize(_nig_cdf_scalar, otypes=[float])


def _scipy_to_wiki(a, b, loc, scale):
    """Convert SciPy (a, b, loc, scale) to Wikipedia (alpha, beta, mu, delta)."""
    return a / scale, b / scale, loc, scale


class _NormalInverseGaussian:
    r"""SciPy-compatible Normal Inverse Gaussian distribution object.

    The single instance :data:`nig.norminvgauss` mirrors the interface of
    ``scipy.stats.norminvgauss``.

    Methods
    -------
    cdf(x, a, b, loc=0, scale=1)
        Cumulative distribution function.
    __call__(a, b, loc=0, scale=1)
        Return a *frozen* distribution with fixed parameters.

    Examples
    --------
    Direct call (mirrors SciPy's unbound method):

    >>> from nig import norminvgauss
    >>> norminvgauss.cdf(2.0, a=4.0, b=-0.8, loc=1.75, scale=2.0)

    Frozen distribution (mirrors SciPy's frozen rv):

    >>> dist = norminvgauss(a=4.0, b=-0.8, loc=1.75, scale=2.0)
    >>> dist.cdf(2.0)
    """

    def cdf(self, x, a, b, loc=0.0, scale=1.0):
        """Cumulative distribution function.

        Parameters
        ----------
        x : float or array-like
            Evaluation point(s).
        a : float
            Tail-heaviness parameter (``a = alpha * delta``, ``a > |b|``).
        b : float
            Asymmetry parameter (``b = beta * delta``).
        loc : float, optional
            Location parameter (default 0).
        scale : float, optional
            Scale parameter (default 1, must be > 0).

        Returns
        -------
        float or ndarray
        """
        alpha, beta, mu, delta = _scipy_to_wiki(a, b, loc, scale)
        # The C++ integration kernel computes 1/t near t=0 (valid IEEE 754,
        # yielding inf which is handled internally). numpy's FPE state flags
        # this transiently; suppress the spurious warning here.
        with np.errstate(divide='ignore', invalid='ignore'):
            result = _nig_cdf_vec(x, alpha, beta, mu, delta)
        return result.item() if result.ndim == 0 else result

    def __call__(self, a, b, loc=0.0, scale=1.0):
        """Freeze the distribution with fixed parameters.

        Returns
        -------
        _FrozenNIG
            Object whose methods accept only ``x``.
        """
        return _FrozenNIG(self, a, b, loc, scale)


class _FrozenNIG:
    """Frozen Normal Inverse Gaussian distribution (fixed parameters).

    Obtain via ``norminvgauss(a, b, loc, scale)``.
    """

    def __init__(self, dist, a, b, loc, scale):
        self._dist = dist
        self.a = a
        self.b = b
        self.loc = loc
        self.scale = scale

    def cdf(self, x):
        """CDF at ``x`` using the frozen parameters."""
        return self._dist.cdf(x, self.a, self.b, self.loc, self.scale)


#: Single instance – use this like ``scipy.stats.norminvgauss``.
norminvgauss = _NormalInverseGaussian()
