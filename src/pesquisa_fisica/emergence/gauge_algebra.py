"""Entropia/IM de regiões em teoria de gauge 1+1D (Schwinger, Gauss eliminada) para diferentes ÁLGEBRAS.

Modos (Casini-Huerta-Rosabal 2014):
  'qubit'  álgebra completa dos qubits da região (inclui linhas de Wilson implícitas entre regiões)
  'inv'    invariante sem centro: operadores que comutam com a carga Q_R de cada região R
           (sem linhas de Wilson entre regiões) → ρ decomposto por setores de Q_R
  'ce'     centro elétrico: 'inv' + fluxos elétricos de fronteira; em 1D aberto o fluxo à esquerda de R é
           a carga total à esquerda, L = Σ_{k<min R} Q_k (variável clássica)
  'carga'  só a distribuição de cargas (diagonal na base σᶻ)

Base: bit n = sítio n (1 = ocupado); Q_n = bit_n − (n mod 2).  Eixo k do tensor (2,)*N ↔ sítio N−1−k.
"""
from __future__ import annotations

import numpy as np


def _bits(idx: np.ndarray, nsites: int) -> np.ndarray:
    return (idx[:, None] >> (nsites - 1 - np.arange(nsites))) & 1


def _charge(bits: np.ndarray, sites) -> np.ndarray:
    par = np.array([s % 2 for s in sites])
    return (bits - par).sum(axis=1) if len(sites) else np.zeros(bits.shape[0], int)


def _H(eigs):
    w = np.asarray(eigs)
    w = w[w > 1e-14]
    return float(-(w * np.log(w)).sum())


def algebra_entropy(full: np.ndarray, N: int, regions, classical_sets=(), mode: str = "qubit") -> float:
    sites = [s for R in regions for s in R]
    rest = [k for k in range(N) if k not in sites]
    T = full.reshape((2,) * N).transpose([N - 1 - s for s in sites] + [N - 1 - s for s in rest])
    X = T.reshape(1 << len(sites), -1)
    rb = _bits(np.arange(X.shape[0]), len(sites))
    if mode == "carga":
        return _H((np.abs(X) ** 2).sum(axis=1))
    if mode == "qubit":
        return _H(np.linalg.eigvalsh(X @ X.conj().T))
    off, mask = 0, np.ones((X.shape[0],) * 2, bool)
    for R in regions:
        q = _charge(rb[:, off:off + len(R)], R)
        mask &= q[:, None] == q[None, :]
        off += len(R)
    if mode == "inv":
        return _H(np.linalg.eigvalsh((X @ X.conj().T) * mask))
    if mode == "ce":
        cb = _bits(np.arange(X.shape[1]), len(rest))
        pos = {s: i for i, s in enumerate(rest)}
        key = np.zeros(X.shape[1], np.int64)
        for C in classical_sets:
            cols = [pos[s] for s in C]
            qc = _charge(cb[:, cols], C) if C else np.zeros(X.shape[1], int)
            key = key * (2 * N + 1) + (qc + N)
        eig = []
        for k in np.unique(key):
            Xc = X[:, key == k]
            eig.append(np.linalg.eigvalsh((Xc @ Xc.conj().T) * mask))
        return _H(np.concatenate(eig))
    raise ValueError(mode)


def algebra_mi(full, N, A, B, mode):
    """I(A:B) na álgebra dada; A à esquerda de B, blocos contíguos."""
    A, B = sorted(A), sorted(B)
    left_A = list(range(0, A[0]))
    left_B = list(range(0, B[0]))
    gap = list(range(A[-1] + 1, B[0]))
    cA = [left_A] if mode == "ce" else []
    cB = [left_B] if mode == "ce" else []
    cAB = [left_A, gap] if mode == "ce" else []
    return (algebra_entropy(full, N, [A], cA, mode) + algebra_entropy(full, N, [B], cB, mode)
            - algebra_entropy(full, N, [A, B], cAB, mode))
