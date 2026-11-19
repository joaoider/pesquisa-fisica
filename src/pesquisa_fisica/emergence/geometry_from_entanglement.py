"""Reconstrução de geometria a partir da informação mútua (Cao, Carroll & Michalakis 2017).

Ideia: graus de liberdade que compartilham muita informação mútua estão "perto".
    d_ij = -ln( I_ij / I_max )
Depois, um MDS clássico (escalonamento multidimensional) mostra quantas dimensões
são necessárias para acomodar essas distâncias.

Nenhuma posição espacial entra no cálculo: só o estado quântico.
"""
from __future__ import annotations

import numpy as np


def emergent_distances(I: np.ndarray, floor: float = 1e-300) -> np.ndarray:
    Imax = I.max()
    if Imax <= 0:
        raise ValueError("informação mútua nula: o estado não define geometria (estado produto)")
    D = -np.log(np.maximum(I, floor) / Imax)
    np.fill_diagonal(D, 0.0)
    return D


def classical_mds(D: np.ndarray):
    """Devolve (autovalores em ordem decrescente, coordenadas) do MDS clássico."""
    n = D.shape[0]
    J = np.eye(n) - np.ones((n, n)) / n
    B = -0.5 * J @ (D**2) @ J
    w, V = np.linalg.eigh(B)
    order = np.argsort(w)[::-1]
    w, V = w[order], V[:, order]
    coords = V * np.sqrt(np.clip(w, 0, None))
    return w, coords


def variance_fraction(eigvals: np.ndarray, k: int) -> float:
    pos = np.clip(eigvals, 0, None)
    return float(pos[:k].sum() / pos.sum())


def spearman(a: np.ndarray, b: np.ndarray) -> float:
    ra = np.argsort(np.argsort(a)).astype(float)
    rb = np.argsort(np.argsort(b)).astype(float)
    return float(np.corrcoef(ra, rb)[0, 1])


def true_distances(shape: tuple[int, ...]) -> np.ndarray:
    pts = np.array(np.unravel_index(np.arange(int(np.prod(shape))), shape)).T.astype(float)
    diff = pts[:, None, :] - pts[None, :, :]
    return np.sqrt((diff**2).sum(-1))


def upper(M: np.ndarray) -> np.ndarray:
    return M[np.triu_indices_from(M, k=1)]
