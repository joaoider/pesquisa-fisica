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


def heat_trace_dirac_squared_S4(s: np.ndarray, nmax: int = 4000) -> np.ndarray:
    """Tr e^{-s D̸²} para o operador de Dirac na S⁴ unitária (spinor de 4 componentes).

    Autovalores de D̸: ±(n + 2), cada sinal com degenerescência 4·C(n+3, 3).
    D̸² = -∇² + R/4 (Lichnerowicz), logo E = -R/4.
    """
    s = np.atleast_1d(np.asarray(s, dtype=float))
    n = np.arange(nmax + 1, dtype=float)
    lam = (n + 2) ** 2
    deg = 2 * 4 * (n + 1) * (n + 2) * (n + 3) / 6
    return np.array([np.sum(deg * np.exp(-si * lam)) for si in s])


def heat_trace_hodge_1forms_S4(s: np.ndarray, lmax: int = 4000) -> np.ndarray:
    """Tr e^{-s Δ₁} do Laplaciano de Hodge-de Rham em 1-formas na S⁴ unitária.

    Δ₁ = -∇² + Ric (Weitzenböck). Espectro = co-exatas + exatas:
      co-exatas: (l+1)(l+2), deg l(l+3)(2l+3)/2, l ≥ 1   (l = 1 são os 10 vetores de Killing)
      exatas (dφ): l(l+3), deg igual à dos escalares, l ≥ 1
    """
    s = np.atleast_1d(np.asarray(s, dtype=float))
    l = np.arange(1, lmax + 1, dtype=float)
    lam_co, deg_co = (l + 1) * (l + 2), l * (l + 3) * (2 * l + 3) / 2
    lam_ex, deg_ex = l * (l + 3), (2 * l + 3) * (l + 1) * (l + 2) / 6
    return np.array([np.sum(deg_co * np.exp(-si * lam_co)) + np.sum(deg_ex * np.exp(-si * lam_ex)) for si in s])


def heat_trace_gauge_vector_S4(s: np.ndarray, lmax: int = 4000) -> np.ndarray:
    """Vetor de gauge sem massa (calibre de Feynman): Tr_1formas - 2 Tr_escalares (fantasmas)."""
    return heat_trace_hodge_1forms_S4(s, lmax) - 2 * heat_trace_sphere(s, 4, 0.0, lmax)


def _fit(trace_fn, dim, n_coeffs, s_min, s_max, n_pts):
    s = np.linspace(s_min, s_max, n_pts)
    y = trace_fn(s) * (4 * np.pi * s) ** (dim / 2) / sphere_volume(dim)
    V = np.vander(s, n_coeffs + 2, increasing=True)
    coeffs, *_ = np.linalg.lstsq(V, y, rcond=None)
    return coeffs[:n_coeffs]


def fit_heat_coefficients(dim: int, xi: float = 0.0, n_coeffs: int = 4,
                          s_min: float = 0.002, s_max: float = 0.05, n_pts: int = 60):
    """Ajusta (4πs)^{d/2} Tr e^{-sD} / Vol = Σ_k b_k s^k e devolve b_0..b_{n-1}.

    Para espaços homogêneos, b_k = tr a_k (densidade constante).
    """
    lmax = int(np.sqrt(60 / s_min)) + 50
    return _fit(lambda s: heat_trace_sphere(s, dim, xi, lmax=lmax), dim, n_coeffs, s_min, s_max, n_pts)


def fit_spin_coefficients(field: str, n_coeffs: int = 3, s_min: float = 0.004, s_max: float = 0.06, n_pts: int = 60):
    """Coeficientes b_k na S⁴ unitária para field ∈ {'dirac', 'hodge1', 'gauge_vector'}."""
    lmax = int(np.sqrt(60 / s_min)) + 50
    fn = {"dirac": heat_trace_dirac_squared_S4,
          "hodge1": heat_trace_hodge_1forms_S4,
          "gauge_vector": heat_trace_gauge_vector_S4}[field]
    return _fit(lambda s: fn(s, lmax), 4, n_coeffs, s_min, s_max, n_pts)
