"""Hamiltoniano de emaranhamento (modular) e primeira lei do emaranhamento (H2 → Einstein).

Caminho de Jacobson (1995, 2015) e Faulkner-Guica-Hartman-Myers-Van Raamsdonk (2014): se
  (i)  o Hamiltoniano de emaranhamento K do vácuo restrito a uma região é LOCAL e GEOMÉTRICO
       (Bisognano-Wichmann: K = ∫ β(x) T₀₀(x) dx, com β(x) = 2π (x)(ℓ−x)/(ℓ v) num intervalo de
       comprimento ℓ: a "temperatura de Unruh" local do vácuo), e
  (ii) vale a primeira lei do emaranhamento δS = δ⟨K⟩ para perturbações pequenas,
então, num contexto holográfico, impor (ii) para todas as regiões equivale às equações de Einstein
linearizadas. Aqui testamos (i) e (ii) numa matéria de rede: cadeia crítica de férmions livres
(CFT com c = 1, v = 2w sin k_F = 2 a meio preenchimento).

Para estados gaussianos: ρ_A ∝ exp(−Σ h_ij c_i† c_j) com h = ln[(1 − C_A)/C_A] (C_A = ⟨c_i† c_j⟩ em A).
O logaritmo é mal condicionado (autovalores exponencialmente próximos de 0 e 1): usamos mpmath.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np


def infinite_chain_correlation(ell: int, kF: float = np.pi / 2, dps: int = 60):
    """C_ij = sin(kF (i−j)) / (π (i−j)), C_ii = kF/π (cadeia infinita, exata), em precisão mp."""
    mp.mp.dps = dps
    kF = mp.mpf(kF)
    C = mp.matrix(ell, ell)
    for i in range(ell):
        for j in range(ell):
            C[i, j] = kF / mp.pi if i == j else mp.sin(kF * (i - j)) / (mp.pi * (i - j))
    return C


def entanglement_hamiltonian(C_mp) -> np.ndarray:
    """h = ln[(1 − C)/C] calculado em precisão estendida; devolve float64."""
    E, Q = mp.eighe(C_mp) if hasattr(mp, "eighe") else mp.eigsy(C_mp)
    n = C_mp.rows
    D = mp.matrix(n, n)
    for k in range(n):
        lam = E[k]
        D[k, k] = mp.log((1 - lam) / lam)
    H = Q * D * Q.T
    return np.array(H.tolist(), dtype=float)


def bw_profile(ell: int, v: float = 2.0) -> np.ndarray:
    """β nas ligações (x + ½), x = 0..ℓ−2, com o intervalo ocupando [−½, ℓ − ½]."""
    x = np.arange(ell - 1) + 1.0      # posição da ligação medida a partir da borda esquerda (−½)
    return 2 * np.pi * x * (ell - x) / (ell * v)


# ---------------------------------------------------------------- primeira lei em cadeia finita

def chain_correlation(N: int, potential=None, w: float = 1.0) -> np.ndarray:
    h = np.zeros((N, N))
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = -w
    if potential is not None:
        h += np.diag(potential)
    e, U = np.linalg.eigh(h)
    occ = U[:, : N // 2]
    return occ @ occ.T


def fermion_entropy(C: np.ndarray) -> float:
    lam = np.clip(np.linalg.eigvalsh(C), 1e-15, 1 - 1e-15)
    return float(-np.sum(lam * np.log(lam) + (1 - lam) * np.log(1 - lam)))


def first_law_check(N: int, ell: int, potential: np.ndarray, eps: float, v: float = 2.0):
    """Devolve (δS, δ⟨K_BW⟩) para a perturbação ε·V do Hamiltoniano, intervalo central de comprimento ℓ.
    K_BW = Σ_ligações β(x+½)·(−(c_x† c_{x+1} + h.c.))  (densidade de energia de ligação, sem termo de μ
    porque o meio preenchimento é simétrico partícula-buraco)."""
    a = (N - ell) // 2
    A = slice(a, a + ell)
    C0 = chain_correlation(N)
    C1 = chain_correlation(N, eps * potential)
    dS = fermion_entropy(C1[A, A]) - fermion_entropy(C0[A, A])
    beta = bw_profile(ell, v)
    bond = lambda C: np.array([-2 * C[a + x, a + x + 1] for x in range(ell - 1)])
    dK = float(np.sum(beta * (bond(C1) - bond(C0))))
    return dS, dK


def chain_correlation_bonds(N: int, bond_mod=None, w: float = 1.0) -> np.ndarray:
    """Cadeia com salto modulado w_x = w (1 + f_x): perturbação da 'métrica' local (velocidade local)."""
    h = np.zeros((N, N))
    f = np.zeros(N - 1) if bond_mod is None else bond_mod
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = -w * (1 + f[i])
    e, U = np.linalg.eigh(h)
    occ = U[:, : N // 2]
    return occ @ occ.T


def first_law_bonds(N: int, ell: int, f: np.ndarray, eps: float, K_exact: np.ndarray | None = None, v: float = 2.0):
    """(δS, δ⟨K_BW⟩, δ⟨K_exato⟩) para a modulação de salto ε·f, intervalo central de comprimento ℓ.
    K_exato: Hamiltoniano de emaranhamento de partícula única do vácuo (ℓ×ℓ), δ⟨K⟩ = Σ h_ij δC_ij."""
    a = (N - ell) // 2
    A = slice(a, a + ell)
    C0 = chain_correlation_bonds(N)
    C1 = chain_correlation_bonds(N, eps * f)
    dS = fermion_entropy(C1[A, A]) - fermion_entropy(C0[A, A])
    beta = bw_profile(ell, v)
    bond = lambda C: np.array([-2 * C[a + x, a + x + 1] for x in range(ell - 1)])
    dK_bw = float(np.sum(beta * (bond(C1) - bond(C0))))
    dK_ex = None if K_exact is None else float(np.sum(K_exact * (C1[A, A] - C0[A, A])))
    return dS, dK_bw, dK_ex
