"""Gauge Z₂ pura em 2+1D via dualidade exata com o Ising de campo transverso na rede DUAL (plaquetes).

H_gauge = −Σ_p B_p − h Σ_ℓ σˣ_ℓ  (bordas lisas: só elos internos ao conjunto de plaquetes)
      ⟷  H_dual = −Σ_p μˣ_p − h Σ_⟨pq⟩ μᶻ_p μᶻ_q
B_p = μˣ_p (plaquete, invariante de gauge) ; σˣ_ℓ = μᶻ_p μᶻ_q (campo elétrico no elo comum).
Um μᶻ_p isolado é uma CORDA elétrica de p até a borda (não local no gauge).
Trabalhamos na base de μˣ (bit 1 ↔ μˣ = −1) e no setor par (Π μˣ = +1). Sítio s = x·Ly + y; bit s = sítio s.
Desconfinado (h pequeno) ↔ paramagneto dual; confinado (h grande) ↔ ordenado dual; h_c ≈ 0,3285 (2D TFIM).
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def bonds(Lx: int, Ly: int):
    s = lambda x, y: x * Ly + y
    out = []
    for x in range(Lx):
        for y in range(Ly):
            if x + 1 < Lx:
                out.append((s(x, y), s(x + 1, y)))
            if y + 1 < Ly:
                out.append((s(x, y), s(x, y + 1)))
    return out


def even_basis(n: int) -> np.ndarray:
    idx = np.arange(1 << n, dtype=np.int64)
    par = np.zeros_like(idx)
    for k in range(n):
        par ^= (idx >> k) & 1
    return idx[par == 0]


def hamiltonian(Lx: int, Ly: int, h: float):
    n = Lx * Ly
    basis = even_basis(n)
    bits = ((basis[:, None] >> np.arange(n)) & 1)
    diag = -(1 - 2 * bits).sum(1).astype(float)
    rows = [np.arange(len(basis))]
    cols = [np.arange(len(basis))]
    vals = [diag]
    for (i, j) in bonds(Lx, Ly):
        tgt = basis ^ ((1 << i) | (1 << j))
        rows.append(np.arange(len(basis)))
        cols.append(np.searchsorted(basis, tgt))
        vals.append(np.full(len(basis), -h))
    H = coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                   shape=(len(basis),) * 2).tocsr()
    return H, basis


def ground_state(Lx: int, Ly: int, h: float):
    H, basis = hamiltonian(Lx, Ly, h)
    e, V = eigsh(H, k=2, which="SA", tol=1e-12)
    o = np.argsort(e)
    full = np.zeros(1 << (Lx * Ly))
    full[basis] = V[:, o[0]]
    return float(e[o[0]]), float(e[o[1]] - e[o[0]]), full


def _H(p):
    p = np.asarray(p, float)
    p = p[p > 1e-15]
    return float(-(p * np.log(p)).sum())


def mi_matrices(full: np.ndarray, n: int):
    """IM entre plaquetes únicos: 'full' (álgebra completa dos spins duais) e 'inv' (só B_p: diagonal em μˣ)."""
    from .ed_quench import rdm
    S1f, S1i = np.zeros(n), np.zeros(n)
    for p in range(n):
        r = rdm(full, n, [p])
        S1f[p] = _H(np.linalg.eigvalsh(r)); S1i[p] = _H(np.real(np.diag(r)))
    If, Ii = np.zeros((n, n)), np.zeros((n, n))
    for p in range(n):
        for q in range(p + 1, n):
            r = rdm(full, n, [p, q])
            If[p, q] = If[q, p] = S1f[p] + S1f[q] - _H(np.linalg.eigvalsh(r))
            Ii[p, q] = Ii[q, p] = S1i[p] + S1i[q] - _H(np.real(np.diag(r)))
    return If, Ii
