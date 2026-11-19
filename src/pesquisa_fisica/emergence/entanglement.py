"""Emaranhamento em estados gaussianos (H2: o espaço emerge do emaranhamento).

Dois modelos de brinquedo de MATÉRIA livre na rede:
  - cadeia de férmions livres (tight-binding) no estado fundamental;
  - rede harmônica (campo escalar discretizado) com massa m, em d = 1 ou 2.

Para estados gaussianos a entropia de qualquer região sai das matrizes de correlação
(Peschel 2003; Casini & Huerta 2009), então dá para calcular redes de centenas de sítios exatamente.
"""
from __future__ import annotations

import numpy as np


# ---------------------------------------------------------------- férmions livres

def fermion_chain_correlations(n_sites: int, periodic: bool = True) -> np.ndarray:
    """C_ij = <c_i^† c_j> no estado fundamental de H = -Σ (c_i^† c_{i+1} + h.c.) a meio preenchimento."""
    H = np.zeros((n_sites, n_sites))
    for i in range(n_sites - 1):
        H[i, i + 1] = H[i + 1, i] = -1.0
    if periodic:
        H[0, -1] = H[-1, 0] = -1.0
    eps, U = np.linalg.eigh(H)
    occ = U[:, eps < 0]
    return occ @ occ.T


def fermion_entropy(C: np.ndarray, region) -> float:
    lam = np.linalg.eigvalsh(C[np.ix_(region, region)])
    lam = np.clip(lam, 1e-14, 1 - 1e-14)
    return float(-np.sum(lam * np.log(lam) + (1 - lam) * np.log(1 - lam)))


# ---------------------------------------------------------------- rede harmônica (escalar livre)

def lattice_laplacian(shape: tuple[int, ...], periodic: bool = False) -> np.ndarray:
    n = int(np.prod(shape))
    idx = np.arange(n).reshape(shape)
    L = np.zeros((n, n))
    for axis in range(len(shape)):
        a = idx
        b = np.roll(idx, -1, axis=axis)
        pairs = zip(a.ravel(), b.ravel())
        if not periodic:
            sl = [slice(None)] * len(shape)
            sl[axis] = slice(0, shape[axis] - 1)
            pairs = zip(a[tuple(sl)].ravel(), b[tuple(sl)].ravel())
        for i, j in pairs:
            L[i, i] += 1
            L[j, j] += 1
            L[i, j] -= 1
            L[j, i] -= 1
    return L


def harmonic_correlations(shape: tuple[int, ...], mass: float, periodic: bool = False):
    """X = <x x> = ½ K^{-1/2}, P = <p p> = ½ K^{1/2}, com K = m² + (-Laplaciano)."""
    K = mass**2 * np.eye(int(np.prod(shape))) + lattice_laplacian(shape, periodic)
    w, V = np.linalg.eigh(K)
    X = 0.5 * (V / np.sqrt(w)) @ V.T
    P = 0.5 * (V * np.sqrt(w)) @ V.T
    return X, P


def _h(nu: np.ndarray) -> np.ndarray:
    nu = np.maximum(nu, 0.5 + 1e-13)
    return (nu + 0.5) * np.log(nu + 0.5) - (nu - 0.5) * np.log(nu - 0.5)


def harmonic_entropy(X: np.ndarray, P: np.ndarray, region) -> float:
    r = np.asarray(region)
    M = X[np.ix_(r, r)] @ P[np.ix_(r, r)]
    nu = np.sqrt(np.clip(np.linalg.eigvals(M).real, 0.25, None))
    return float(np.sum(_h(nu)))


def harmonic_mutual_information_matrix(X: np.ndarray, P: np.ndarray) -> np.ndarray:
    """I(i:j) = S_i + S_j - S_ij para todos os pares de sítios (estado gaussiano, exato)."""
    n = X.shape[0]
    s1 = np.array([float(_h(np.sqrt(np.array([max(X[i, i] * P[i, i], 0.25)])))[0]) for i in range(n)])
    I = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            I[i, j] = I[j, i] = max(s1[i] + s1[j] - harmonic_entropy(X, P, [i, j]), 0.0)
    return I
