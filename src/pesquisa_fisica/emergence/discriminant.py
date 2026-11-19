"""Teste discriminante: o que organiza a geometria emergente é só o GAP, ou algo específico da FORÇA?

Hipótese nula (literatura: gap ⇒ comprimento de correlação ⇒ informação mútua local; "entanglement não
distingue confinamento de gap"): com o MESMO gap, o perfil de informação mútua I(d) é o mesmo, venha o gap
do acoplamento de gauge J (confinamento) ou de uma massa de férmion sem gauge.

Protocolo:
  1. caso de gauge (m_q, J) → gap Δ = E₁ − E₀ (setor físico), perfil I(d) na partição fermiônica;
  2. caso livre (J = 0) com massa staggered m' ajustada por bissecção para o MESMO Δ;
  3. comparar I_gauge(d)/I_livre(d).
Controle adicional: igualar a taxa de decaimento assintótica (comprimento de correlação) em vez do gap.

Mecanismo esperado se a nula cair (bosonização): na teoria confinada as excitações com gap são MÉSONS
(compostos bosônicos); na livre, o portador mais leve é um férmion com metade do gap do par. Como I ~ (correlação)²,
no contínuo o decaimento confinado seria ~2× mais rápido no mesmo gap.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh

from . import schwinger as sw
from .su2_lattice import SU2Chain
from .su2_lattice import site_mutual_information as fermionic_mi


def profile_from_matrix(I: np.ndarray) -> np.ndarray:
    n = len(I)
    return np.array([np.mean([I[i, i + d] for i in range(n - d)]) for d in range(1, n)])


# ---------------------------------------------------------------- U(1): Schwinger

def schwinger_gap_and_profile(N: int, m: float, J: float):
    H, basis = sw.hamiltonian(N, w=1.0, m=m, J=J)
    w, V = eigsh(H, k=2, which="SA", tol=1e-12)
    o = np.argsort(w)
    v = V[:, o[0]]
    I = fermionic_mi(basis, v, N, modes_per_site=2)            # células = pares (elétron, pósitron)
    return float(w[o[1]] - w[o[0]]), profile_from_matrix(I)


def schwinger_gap(N: int, m: float, J: float) -> float:
    H, _ = sw.hamiltonian(N, w=1.0, m=m, J=J)
    w = np.sort(eigsh(H, k=2, which="SA", tol=1e-12)[0])
    return float(w[1] - w[0])


# ---------------------------------------------------------------- SU(2)

def su2_gap_and_profile(N: int, m: float, J: float):
    ch = SU2Chain(N, w=1.0, m=m, J=J)
    op = LinearOperator((ch.dim, ch.dim), matvec=ch.matvec, dtype=float)
    w, V = eigsh(op, k=2, which="SA", tol=1e-11, v0=np.random.default_rng(0).standard_normal(ch.dim))
    o = np.argsort(w)
    v = V[:, o[0]]
    return float(w[o[1]] - w[o[0]]), profile_from_matrix(fermionic_mi(ch.basis, v, 2 * N)), ch.total_color_casimir(v)


# ---------------------------------------------------------------- pareamento

def match_mass_to_gap(gap_fn, target: float, lo: float = 0.0, hi: float = 8.0, iters: int = 22) -> float:
    """Bissecção em m' (com J = 0) até gap_fn(m') = target (gap cresce com m')."""
    for _ in range(iters):
        mid = (lo + hi) / 2
        if gap_fn(mid) < target:
            lo = mid
        else:
            hi = mid
    return (lo + hi) / 2


def decay_rate(profile: np.ndarray, d_from: int = 3, d_to: int = 5) -> float:
    d = np.arange(d_from, d_to + 1)
    return float(-np.polyfit(d, np.log(profile[d_from - 1: d_to]), 1)[0])
