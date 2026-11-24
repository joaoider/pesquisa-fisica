"""Gauge Z_N pura 2+1D, rede aberta Px×Py com todos os elos. Dual exato em alturas n_p ∈ Z_N nos plaquetes.

H = (h/2) Σ_ℓ s(E_ℓ)² − K Σ_p (U_p + U_p†)/2,   s(E) = representante simétrico de E mod N.
Elo orientado v1→v2 (+x̂ ou +ŷ): E_ℓ = n_esq − n_dir (mod N), plaquete ausente = 0 ⇒ Gauss automática.
U_p (magnético) desloca n_p → n_p + 1.  N = 2 reproduz o Z₂ com h = 4·h_Z₂.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh

from .z2_links import Z2Lattice


class ZNLattice(Z2Lattice):
    def __init__(self, Px: int, Py: int, N: int):
        super().__init__(Px, Py)
        self.N = N
        pid = lambda x, y: x * Py + y if (0 <= x < Px and 0 <= y < Py) else -1
        self.lr = []
        for (v1, v2, _) in self.links:
            (x, y), (x2, y2) = v1, v2
            if y2 == y:                       # horizontal +x: esquerda = acima, direita = abaixo
                self.lr.append((pid(x, y), pid(x, y - 1)))
            else:                             # vertical +y: esquerda = −x, direita = +x
                self.lr.append((pid(x - 1, y), pid(x, y)))
        self.dim = N ** self.P

    def digits(self) -> np.ndarray:
        c = np.arange(self.dim, dtype=np.int64)
        return np.stack([(c // self.N ** p) % self.N for p in range(self.P)], axis=1)

    def link_values(self, dig: np.ndarray) -> np.ndarray:
        n = lambda p: dig[:, p] if p >= 0 else 0
        return np.stack([(n(l) - n(r)) % self.N for (l, r) in self.lr], axis=1)

    def sym(self, E):
        return np.where(E <= self.N // 2, E, E - self.N)

    def hamiltonian(self, h: float, K: float = 1.0):
        dig = self.digits()
        E = self.link_values(dig)
        diag = 0.5 * h * (self.sym(E) ** 2).sum(1).astype(float)
        c = np.arange(self.dim, dtype=np.int64)
        rows, cols, vals = [c], [c], [diag]
        for p in range(self.P):
            step = self.N ** p
            up = np.where(dig[:, p] < self.N - 1, c + step, c - (self.N - 1) * step)
            dn = np.where(dig[:, p] > 0, c - step, c + (self.N - 1) * step)
            for tgt in (up, dn):
                rows.append(c); cols.append(tgt); vals.append(np.full(self.dim, -K / 2))
        return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))),
                          shape=(self.dim,) * 2).tocsr()

    def ground_state(self, h: float, K: float = 1.0):
        e, V = eigsh(self.hamiltonian(h, K), k=2, which="SA", tol=1e-10)
        o = np.argsort(e)
        psi = V[:, o[0]].real
        return float(e[o[0]]), float(e[o[1]] - e[o[0]]), psi / np.linalg.norm(psi)


def _H(p):
    p = np.asarray(p, float)
    p = p[p > 1e-15]
    return float(-(p * np.log(p)).sum())


def region_entropies_zn(lat: ZNLattice, psi, A, E, max_dim: int = 6000):
    N = lat.N
    Aset = set(A)
    Ac = [k for k in range(lat.L) if k not in Aset]
    enc = lambda cols: (E[:, cols] * (N ** np.arange(len(cols)))).sum(1) if cols else np.zeros(len(psi), np.int64)
    a, b = enc(A), enc(Ac)
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    if min(len(ua), len(ub)) > max_dim:
        return None
    M = coo_matrix((psi, (ia, ib)), shape=(len(ua), len(ub))).tocsr()
    G = (M @ M.T) if len(ua) <= len(ub) else (M.T @ M)
    S_ext = _H(np.linalg.eigvalsh(G.toarray()))
    cuts, inc = lat.cut_vertices(A)
    f = np.zeros(len(psi), np.int64)
    for j, v in enumerate(cuts):
        flux = np.zeros(len(psi), np.int64)
        for k in inc[v]:
            if k in Aset:
                v1, v2, _ = lat.links[k]
                flux += E[:, k] if v1 == v else -E[:, k]
        f += (flux % N) * N ** j
    _, inv = np.unique(f, return_inverse=True)
    H_edge = _H(np.bincount(inv, weights=psi ** 2))
    return {"S_ext": S_ext, "H_edge": H_edge, "S_dist": S_ext - H_edge, "n_cut": len(cuts), "nA": len(A)}
