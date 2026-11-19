"""Modelo de Kitaev no favo de mel: teoria de gauge Z2 + matéria (Majorana), exata em redes grandes.

H = -Σ_α J_α Σ_{<jk>_α} σ_j^α σ_k^α.  Com σ^α = i b^α c (Kitaev 2006):
    H = (i/4) Σ_{jk} A_jk c_j c_k,   A_jk = 2 J_α u_jk  (j ∈ A, k ∈ B),  u_jk = i b_j^α b_k^α = ±1
  - u_jk: campo de gauge Z2 nas ligações (estático; fluxo W_p = Π u em cada hexágono);
  - c_j: matéria fermiônica de Majorana, livre em cada configuração de u;
  - transformação de gauge no sítio j: u → -u nas 3 ligações de j  ⇔  c_j → -c_j.
Fase sem vórtices (u = +1 uniforme) contém o estado fundamental (teorema de Lieb).
Fases: gapless B (|J_x|, |J_y|, |J_z| obedecem desigualdade triangular), gapped A (J_z > J_x + J_y).

O emaranhamento total de uma região é S = S_matéria + S_gauge, com S_gauge = (n_ligações cortadas - 1) ln 2
(Yao & Qi 2010): a parte de gauge só conta fronteira. A geometria vem da matéria.
"""
from __future__ import annotations

from dataclasses import dataclass, field
from functools import cached_property

import numpy as np

A1 = np.array([1.0, 0.0])
A2 = np.array([0.5, np.sqrt(3) / 2])


@dataclass
class Honeycomb:
    L: int
    J: tuple = (1.0, 1.0, 1.0)
    u: np.ndarray | None = field(default=None, repr=False)   # u por ligação (3·L² ligações)

    @property
    def n_cells(self) -> int:
        return self.L * self.L

    def cell(self, n1: int, n2: int) -> int:
        return (n1 % self.L) + self.L * (n2 % self.L)

    @cached_property
    def bonds(self) -> list[tuple[int, int, int]]:
        """(sítio A, sítio B, tipo α). Sítio A da célula r = 2r, sítio B = 2r+1.
        x: A(r)–B(r); y: A(r)–B(r - a1); z: A(r)–B(r - a2)."""
        out = []
        for n2 in range(self.L):
            for n1 in range(self.L):
                r = self.cell(n1, n2)
                out.append((2 * r, 2 * r + 1, 0))
                out.append((2 * r, 2 * self.cell(n1 - 1, n2) + 1, 1))
                out.append((2 * r, 2 * self.cell(n1, n2 - 1) + 1, 2))
        return out

    def gauge_field(self) -> np.ndarray:
        return np.ones(len(self.bonds)) if self.u is None else self.u

    def A_matrix(self) -> np.ndarray:
        n = 2 * self.n_cells
        A = np.zeros((n, n))
        for (j, k, a), u in zip(self.bonds, self.gauge_field()):
            A[j, k] += 2 * self.J[a] * u
            A[k, j] -= 2 * self.J[a] * u
        return A

    def ground_state_covariance(self) -> tuple[float, np.ndarray]:
        """Γ_jk = (i/2)<[c_j, c_k]> real antissimétrica do estado fundamental da matéria; devolve (E, Γ)."""
        A = self.A_matrix()
        w, V = np.linalg.eigh(1j * A)
        sgn = (V * np.sign(w)) @ V.conj().T
        best = None
        for s in (+1, -1):
            G = np.real(-1j * s * sgn)
            E = -0.25 * np.sum(A * G)
            if best is None or E < best[0]:
                best = (E, G)
        return best

    @cached_property
    def cell_positions(self) -> np.ndarray:
        n1, n2 = np.meshgrid(np.arange(self.L), np.arange(self.L), indexing="xy")
        return n1.ravel()[:, None] * A1 + n2.ravel()[:, None] * A2

    def torus_cell_distances(self) -> np.ndarray:
        L = self.L
        n1, n2 = np.meshgrid(np.arange(L), np.arange(L), indexing="xy")
        n1, n2 = n1.ravel(), n2.ravel()
        d1 = (n1[:, None] - n1[None, :]) % L
        d2 = (n2[:, None] - n2[None, :]) % L
        best = np.full(d1.shape, np.inf)
        for s1 in (-L, 0, L):
            for s2 in (-L, 0, L):
                v = (d1 + s1)[..., None] * A1 + (d2 + s2)[..., None] * A2
                best = np.minimum(best, np.sqrt((v**2).sum(-1)))
        return best


def gauge_transform(lat: Honeycomb, sites) -> np.ndarray:
    """u' após transformação de gauge D_j nos sítios dados (inverte u nas 3 ligações de cada j)."""
    u = lat.gauge_field().copy()
    sites = set(sites)
    for b, (j, k, _) in enumerate(lat.bonds):
        if j in sites:
            u[b] *= -1
        if k in sites:
            u[b] *= -1
    return u


# ---------------------------------------------------------------- entropias de Majorana

def _h2(p):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -(p * np.log(p) + (1 - p) * np.log(1 - p))


def majorana_entropy(G: np.ndarray, modes) -> float:
    m = np.asarray(modes)
    nu = np.abs(np.linalg.eigvals(G[np.ix_(m, m)]).imag)
    nu = np.sort(nu)[::2]  # autovalores ±iν aparecem em pares
    return float(np.sum(_h2((1 + nu) / 2)))


def cell_mutual_information(G: np.ndarray, n_cells: int) -> np.ndarray:
    """I(r:r') entre células (2 Majoranas cada = 1 férmion complexo), vetorizado."""
    s1 = np.array([majorana_entropy(G, [2 * r, 2 * r + 1]) for r in range(n_cells)])
    iu, ju = np.triu_indices(n_cells, 1)
    modes = np.stack([2 * iu, 2 * iu + 1, 2 * ju, 2 * ju + 1], axis=1)
    blocks = G[modes[:, :, None], modes[:, None, :]]
    nu = np.abs(np.linalg.eigvals(blocks).imag)
    nu = np.sort(nu, axis=1)[:, ::2]
    s2 = _h2((1 + nu) / 2).sum(1)
    I = np.zeros((n_cells, n_cells))
    I[iu, ju] = np.maximum(s1[iu] + s1[ju] - s2, 0.0)
    return I + I.T
