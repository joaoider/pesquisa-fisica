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


# ---------------------------------------------------------------- estimadores que só usam a ORDEM da informação mútua
# (não dependem de I decair exponencialmente ou como lei de potência: servem para fases gapped e críticas)

def knn_graph(S: np.ndarray, k: int) -> np.ndarray:
    """Grafo simétrico ligando cada grau de liberdade aos k parceiros de maior S (ex.: informação mútua)."""
    n = len(S)
    S = S.astype(float).copy()
    np.fill_diagonal(S, -np.inf)
    nb = np.argsort(S, axis=1)[:, ::-1][:, :k]
    A = np.zeros((n, n), dtype=bool)
    A[np.repeat(np.arange(n), k), nb.ravel()] = True
    return A | A.T


def graph_distances(A: np.ndarray) -> np.ndarray:
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import shortest_path
    return shortest_path(csr_matrix(A.astype(float)), unweighted=True)


def ball_growth_dimension(D: np.ndarray, radii=(2, 3, 4, 5)) -> float:
    """Inclinação de ln N(r) × ln r, com N(r) = nº médio de nós a distância de grafo ≤ r."""
    N = [(D <= r).sum(1).mean() for r in radii]
    return float(np.polyfit(np.log(radii), np.log(N), 1)[0])
