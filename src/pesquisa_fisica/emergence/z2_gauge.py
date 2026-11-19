"""Teoria de gauge Z2 com matéria na rede (H2 · primeiro passo rumo às "3 forças").

Modelo de gauge-Higgs de Z2 (Fradkin & Shenker 1979) no calibre unitário, em que os campos de
matéria são absorvidos e os qubits das ligações já são graus de liberdade invariantes de gauge:

    H = -Δ Σ_v A_v  -  J Σ_p B_p  -  h Σ_l X_l  -  λ Σ_l Z_l

    A_v = Π_{l ∋ v} X_l    paridade da carga de matéria no vértice v (Δ = massa da matéria)
    B_p = Π_{l ∈ p} Z_l    fluxo magnético na plaqueta p (J = acoplamento de plaqueta)
    X_l                    campo elétrico na ligação (h = tensão elétrica → confinamento)
    Z_l                    salto da matéria τ^z σ^z τ^z ao longo da ligação (λ → fase de Higgs)

Limites: h = λ = 0 é o código tórico (fase desconfinada, ordem topológica);
h grande → confinamento; λ grande → Higgs. Rede quadrada L×L periódica, 2L² ligações.
"""
from __future__ import annotations

from dataclasses import dataclass
from functools import cached_property

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


@dataclass
class SquareTorus:
    L: int

    @property
    def n_links(self) -> int:
        return 2 * self.L * self.L

    def link(self, x: int, y: int, direction: int) -> int:
        """direction 0 = horizontal (x→x+1), 1 = vertical (y→y+1)."""
        L = self.L
        return 2 * ((x % L) + L * (y % L)) + direction

    @cached_property
    def stars(self) -> list[list[int]]:
        L = self.L
        return [[self.link(x, y, 0), self.link(x - 1, y, 0), self.link(x, y, 1), self.link(x, y - 1, 1)]
                for y in range(L) for x in range(L)]

    @cached_property
    def plaquettes(self) -> list[list[int]]:
        L = self.L
        return [[self.link(x, y, 0), self.link(x + 1, y, 1), self.link(x, y + 1, 0), self.link(x, y, 1)]
                for y in range(L) for x in range(L)]

    @cached_property
    def midpoints(self) -> np.ndarray:
        pts = np.zeros((self.n_links, 2))
        for y in range(self.L):
            for x in range(self.L):
                pts[self.link(x, y, 0)] = (x + 0.5, y)
                pts[self.link(x, y, 1)] = (x, y + 0.5)
        return pts

    def torus_distances(self) -> np.ndarray:
        d = np.abs(self.midpoints[:, None, :] - self.midpoints[None, :, :])
        d = np.minimum(d, self.L - d)
        return np.sqrt((d**2).sum(-1))


def _masks(groups: list[list[int]]) -> np.ndarray:
    return np.array([sum(1 << q for q in g) for g in groups], dtype=np.int64)


def ground_state(lat: SquareTorus, Delta=1.0, J=1.0, h=0.0, lam=0.0, seed: int = 0) -> tuple[float, np.ndarray]:
    n = lat.n_links
    dim = 1 << n
    idx = np.arange(dim, dtype=np.int64)
    bits = ((idx[:, None] >> np.arange(n)) & 1).astype(np.int8)           # z_l ∈ {0,1}; Z = 1 - 2z
    zsign = 1 - 2 * bits
    diag = np.zeros(dim)
    for m in _masks(lat.plaquettes):
        parity = (np.bitwise_count(idx & m) & 1).astype(np.int64)
        diag -= J * (1 - 2 * parity)
    diag -= lam * zsign.sum(1)
    del bits, zsign
    star_masks = _masks(lat.stars)
    link_masks = np.array([1 << l for l in range(n)], dtype=np.int64)

    def matvec(v):
        v = np.asarray(v).ravel()
        out = diag * v
        for m in star_masks:
            out -= Delta * v[idx ^ m]
        if h:
            for m in link_masks:
                out -= h * v[idx ^ m]
        return out

    op = LinearOperator((dim, dim), matvec=matvec, dtype=float)
    v0 = np.random.default_rng(seed).standard_normal(dim)
    w, V = eigsh(op, k=1, which="SA", v0=v0, tol=1e-10)
    return float(w[0]), V[:, 0] / np.linalg.norm(V[:, 0])


# ---------------------------------------------------------------- matrizes densidade reduzidas

def reduced_density_matrix(psi: np.ndarray, n: int, keep: list[int]) -> np.ndarray:
    """ρ_A para os qubits em `keep` (bit l do índice = qubit l)."""
    t = psi.reshape([2] * n)                     # eixo 0 = bit n-1, ..., eixo n-1 = bit 0
    axes = [n - 1 - q for q in keep]
    rest = [a for a in range(n) if a not in axes]
    t = np.transpose(t, axes + rest).reshape(1 << len(keep), -1)
    return t @ t.conj().T


def von_neumann(rho: np.ndarray) -> float:
    p = np.linalg.eigvalsh(rho)
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log(p)))


def mutual_information_matrix(psi: np.ndarray, n: int) -> np.ndarray:
    s1 = [von_neumann(reduced_density_matrix(psi, n, [i])) for i in range(n)]
    I = np.zeros((n, n))
    for i in range(n):
        for j in range(i + 1, n):
            I[i, j] = I[j, i] = max(s1[i] + s1[j] - von_neumann(reduced_density_matrix(psi, n, [i, j])), 0.0)
    return I


def expectation_product_x(psi: np.ndarray, n: int, qubits: list[int]) -> float:
    m = sum(1 << q for q in qubits)
    idx = np.arange(1 << n, dtype=np.int64)
    return float(np.dot(psi, psi[idx ^ m]))
