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


# ---------------------------------------------------------------- versão MPS (TeNPy): duas células de 2 sítios

def _occ_bits_tenpy(n_sites: int) -> np.ndarray:
    """Bits de ocupação para índices de ρ vindos de block_pair_rdms (TeNPy SpinHalfSite: p=1 ↔ 'up' = ocupado)."""
    idx = np.arange(1 << n_sites)
    return (idx[:, None] >> (n_sites - 1 - np.arange(n_sites))) & 1


def two_cell_mi(rho16: np.ndarray, cell_starts=(0, 0)) -> dict:
    """IM entre duas células de 2 sítios (ρ sobre a1,a2,b1,b2) nas álgebras qubit, inv e carga.
    cell_starts dá a paridade do 1º sítio de cada célula (Q_k = bit_k − (k mod 2))."""
    from .quench import _ptrace

    b = _occ_bits_tenpy(4)
    par = np.array([cell_starts[0] % 2, (cell_starts[0] + 1) % 2, cell_starts[1] % 2, (cell_starts[1] + 1) % 2])
    q = b - par
    QA, QB = q[:, :2].sum(1), q[:, 2:].sum(1)
    b2 = _occ_bits_tenpy(2)

    def S_masked(r, labels_list):
        m = np.ones(r.shape, bool)
        for L in labels_list:
            m &= L[:, None] == L[None, :]
        return _H(np.linalg.eigvalsh((r + r.conj().T) / 2 * m))

    rA, rB = _ptrace(rho16, [0, 1], 4), _ptrace(rho16, [2, 3], 4)
    qA1 = (b2 - par[:2]).sum(1)
    qB1 = (b2 - par[2:]).sum(1)
    out = {}
    out["qubit"] = S_masked(rA, []) + S_masked(rB, []) - S_masked(rho16, [])
    out["inv"] = S_masked(rA, [qA1]) + S_masked(rB, [qB1]) - S_masked(rho16, [QA, QB])
    diag = lambda r: _H(np.real(np.diag(r)))
    out["carga"] = diag(rA) + diag(rB) - diag(rho16)
    return out
