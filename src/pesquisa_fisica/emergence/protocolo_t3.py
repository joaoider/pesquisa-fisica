"""T3: protocolo experimental realista para a I₃ de Rényi-2 na álgebra invariante (Schwinger digital).

- Trotter Z·(ímpares)·(pares) (ordem validada contra Martinez 2016, v0.30): por passo aplica pares, ímpares, diagonal.
- Desfase: após cada passo, cada qubit sofre Z com probabilidade p (trajetórias).
- Medida: sombras clássicas de 1 disparo (Haar local), estimador em lotes; álgebra invariante = desfasagem por setores
  de carga de cada região (operação linear, aplicada no pós-processamento).
Convenção de bits: bit n = sítio n (1 = ocupado); Q_n = bit_n − (n mod 2).
"""
from __future__ import annotations

import numpy as np

from .ed_quench import bare_vacuum_index, rdm


def diag_energies(N: int, m: float, J: float = 0.0, V: float = 0.0, w_unused=None) -> np.ndarray:
    idx = np.arange(1 << N, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(N)) & 1).astype(float)
    sz = 2 * bits - 1
    stag = (-1.0) ** np.arange(N)
    E = 0.5 * m * (sz * stag).sum(1)
    if J:
        Q = (sz + stag) / 2
        E += J * (np.cumsum(Q, axis=1)[:, :-1] ** 2).sum(1)
    if V:
        E += V * (bits[:, :-1] * bits[:, 1:]).sum(1)
    return E


def _bond_gate(psi, n, c, s, idx):
    mask = (1 << n) | (1 << (n + 1))
    diff = ((idx >> n) & 1) != ((idx >> (n + 1)) & 1)
    partner = idx ^ mask
    out = psi.copy()
    out[diff] = c * psi[diff] - 1j * s * psi[partner[diff]]
    return out


def trotter_trajectory(N: int, Ediag: np.ndarray, dt: float, steps: int, p: float, rng, w: float = 1.0,
                       record=None):
    idx = np.arange(1 << N, dtype=np.int64)
    psi = np.zeros(1 << N, complex)
    psi[bare_vacuum_index(N)] = 1.0
    c, s = np.cos(w * dt), np.sin(w * dt)
    phase = np.exp(-1j * dt * Ediag)
    bitsq = [((idx >> q) & 1) for q in range(N)]
    out = {}
    for k in range(1, steps + 1):
        for n in range(0, N - 1, 2):
            psi = _bond_gate(psi, n, c, s, idx)
        for n in range(1, N - 1, 2):
            psi = _bond_gate(psi, n, c, s, idx)
        psi = psi * phase
        if p > 0:
            for q in np.nonzero(rng.random(N) < p)[0]:
                psi = psi * (1 - 2 * bitsq[q])
        if record is not None and k in record:
            out[k] = psi.copy()
    return psi if record is None else out


def _ptrace_small(rho, keep, n):
    t = rho.reshape([2] * (2 * n))
    drop = [i for i in range(n) if i not in keep]
    for k, i in enumerate(sorted(drop, reverse=True)):
        m = n - k
        t = np.trace(t, axis1=i, axis2=i + m)
    d = 2 ** len(keep)
    return t.reshape(d, d)


def _labels(sites_sub, groups):
    """Rótulos de setor (carga de cada grupo) para a base do sub-sistema com sítios 'sites_sub' (MSB primeiro)."""
    n = len(sites_sub)
    idx = np.arange(1 << n)
    b = (idx[:, None] >> (n - 1 - np.arange(n))) & 1
    lab = np.zeros(1 << n, np.int64)
    for g in groups:
        cols = [sites_sub.index(s) for s in g]
        q = (b[:, cols] - np.array([s % 2 for s in g])).sum(1)
        lab = lab * 32 + (q + 8)
    return lab


SUBSETS = {"A": ["A"], "B": ["B"], "C": ["C"], "AB": ["A", "B"], "BC": ["B", "C"], "AC": ["A", "C"], "ABC": ["A", "B", "C"]}


def I3_renyi2_inv_from_rhoABC(rho, regions: dict, purity_fn=None):
    """I₃ de Rényi-2 na álgebra invariante a partir de ρ sobre (A,B,C) (ordem dos sítios = A+B+C)."""
    allsites = regions["A"] + regions["B"] + regions["C"]
    n = len(allsites)
    S = {}
    for key, regs in SUBSETS.items():
        sites = [s for r in regs for s in regions[r]]
        keep = [allsites.index(s) for s in sites]
        lab = _labels(sites, [regions[r] for r in regs])
        mask = lab[:, None] == lab[None, :]
        if purity_fn is None:
            r = _ptrace_small(rho, keep, n) if len(keep) < n else rho
            pur = float(np.real(np.trace((r * mask) @ (r * mask))))
        else:
            pur = purity_fn(keep, mask)
        S[key] = -np.log(pur)
    return S["A"] + S["B"] + S["C"] - S["AB"] - S["BC"] - S["AC"] + S["ABC"]


def rho_regions(full, N, regions):
    return rdm(full, N, regions["A"] + regions["B"] + regions["C"])


# ---------------------------------------------------------------- sombras clássicas em lotes

def _haar2(rng, k):
    z = (rng.normal(size=(k, 2, 2)) + 1j * rng.normal(size=(k, 2, 2))) / np.sqrt(2)
    q, r = np.linalg.qr(z)
    d = np.diagonal(r, axis1=1, axis2=2)
    return q * (d / np.abs(d))[:, None, :]


def shadow_batches(rho, n_qubits: int, M: int, n_batches: int, rng):
    """Médias de snapshots ρ̂ = ⊗(3 u†|s⟩⟨s|u − 1) por lote, a partir de M disparos (um por unitária)."""
    d = 1 << n_qubits
    sums = np.zeros((n_batches, d, d), complex)
    per = M // n_batches
    I2 = np.eye(2)
    for b in range(n_batches):
        for _ in range(per):
            us = _haar2(rng, n_qubits)
            U = us[0]
            for u in us[1:]:
                U = np.kron(U, u)
            pr = np.clip(np.real(np.einsum("ij,jk,ik->i", U, rho, U.conj())), 0, None)
            s = rng.choice(d, p=pr / pr.sum())
            snap = None
            for q in range(n_qubits):
                bit = (s >> (n_qubits - 1 - q)) & 1
                v = us[q].conj()[bit]                     # ⟨s|u  → u†|s⟩ conj
                proj = np.outer(v, v.conj())              # |u†s⟩⟨u†s| = u†|s⟩⟨s|u
                f = 3 * proj - I2
                snap = f if snap is None else np.kron(snap, f)
            sums[b] += snap
        sums[b] /= per
    return sums


def batch_purity_fn(batches, n):
    def fn(keep, mask):
        R = [(_ptrace_small(B, keep, n) if len(keep) < n else B) for B in batches]
        vals = [np.real(np.trace((R[a] * mask) @ R[b])) for a in range(len(R)) for b in range(len(R)) if a != b]
        return max(float(np.mean(vals)), 1e-6)
    return fn
