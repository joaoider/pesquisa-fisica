"""Coeficientes de heat kernel (Seeley-DeWitt) para operadores do tipo Laplace.

Para D = -(g^{ab} ∇_a ∇_b + E) agindo em campos com n componentes, a expansão é

    Tr e^{-sD} ~ (4π s)^{-d/2} Σ_k s^k ∫ √g tr a_k(x)

com (Vassilevich, Phys. Rept. 388 (2003) 279, eq. 4.26-4.28, sem bordo):

    a_0 = n
    a_1 = (1/6) tr(6E + R)
    a_2 = (1/360) tr(60 □E + 60 R E + 180 E² + 12 □R + 5 R² - 2 Ric² + 2 Riem² + 30 Ω_ab Ω^ab)

Para um escalar real com D = -□ + ξR temos E = -ξR, n = 1, Ω = 0.
"""
from __future__ import annotations

import numpy as np
import sympy as sp

from .geometry import Metric


def scalar_a1(metric: Metric, xi=0):
    """Densidade a_1 para escalar real: (1/6 - ξ) R."""
    return sp.simplify((sp.Rational(1, 6) - xi) * metric.ricci_scalar)


def scalar_a2(metric: Metric, xi=0, include_total_derivatives: bool = True):
    """Densidade a_2 para escalar real com acoplamento ξ."""
    R = metric.ricci_scalar
    E = -xi * R
    expr = (180 * E**2 + 60 * R * E + 5 * R**2 - 2 * metric.ricci_squared + 2 * metric.riemann_squared)
    if include_total_derivatives:
        expr += 60 * metric.box(E) + 12 * metric.box(R)
    return sp.simplify(expr / 360)


# ---------------------------------------------------------------- traços espectrais exatos

def heat_trace_sphere(s: np.ndarray, dim: int, xi: float = 0.0, lmax: int = 4000) -> np.ndarray:
    """Tr e^{-s(-Δ + ξR)} na esfera unitária S^dim, somando o espectro exato.

    Autovalores do Laplaciano em S^d: l(l+d-1); degenerescência
    (2l+d-1)(l+d-2)! / (l! (d-1)!).  R = d(d-1).
    """
    from math import comb

    s = np.atleast_1d(np.asarray(s, dtype=float))
    l = np.arange(lmax + 1)
    lam = l * (l + dim - 1) + xi * dim * (dim - 1)
    if dim == 2:
        deg = 2 * l + 1
    else:
        deg = np.array([(2 * k + dim - 1) * comb(k + dim - 2, k) / (dim - 1) for k in l], dtype=float)
    return np.array([np.sum(deg * np.exp(-si * lam)) for si in s])


def sphere_volume(dim: int) -> float:
    from math import gamma, pi
    return 2 * pi ** ((dim + 1) / 2) / gamma((dim + 1) / 2)


def fit_heat_coefficients(dim: int, xi: float = 0.0, n_coeffs: int = 4,
                          s_min: float = 0.002, s_max: float = 0.05, n_pts: int = 60):
    """Ajusta (4πs)^{d/2} Tr e^{-sD} / Vol = Σ_k b_k s^k e devolve b_0..b_{n-1}.

    Para espaços homogêneos, b_k = a_k (densidade constante).
    """
    s = np.linspace(s_min, s_max, n_pts)
    lmax = int(np.sqrt(60 / s_min)) + 50
    y = heat_trace_sphere(s, dim, xi, lmax=lmax) * (4 * np.pi * s) ** (dim / 2) / sphere_volume(dim)
    V = np.vander(s, n_coeffs + 2, increasing=True)
    coeffs, *_ = np.linalg.lstsq(V, y, rcond=None)
    return coeffs[:n_coeffs]
