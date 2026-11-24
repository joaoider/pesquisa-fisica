"""Gauge Z₂ pura 2+1D em rede ABERTA de Px×Py plaquetes, com TODOS os elos (Gauss em todos os vértices).

Dual exato: σˣ_ℓ = Π μᶻ dos plaquetes adjacentes a ℓ (2 nos elos internos, 1 nos de borda);
H = −Σ_p B_p − h Σ_ℓ σˣ_ℓ  ⟷  H_dual = −Σ_p μˣ_p − h Σ_ℓ Π_{p∋ℓ} μᶻ_p  (sem simetria Z₂: campo de borda).
Configurações μᶻ (2^{PxPy}) ↔ configurações de laços elétricos fechados (base σˣ dos elos) — bijeção.
Hilbert ESTENDIDO de uma região de elos A: ρ_A na base elétrica. Centro elétrico: f_v = Π_{ℓ∈A, ℓ∋v} σˣ_ℓ nos
vértices de corte v ⇒ ρ_A é bloco-diagonal em {f_v} e S_ext = H(f) + Σ p_f S(ρ_f) = H_borda + S_dest.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


class Z2Lattice:
    def __init__(self, Px: int, Py: int):
        self.Px, self.Py = Px, Py
        self.P = Px * Py
        pid = lambda x, y: x * Py + y
        links = []                                   # (vértice1, vértice2, [plaquetes])
        for x in range(Px):
            for y in range(Py + 1):                  # horizontais
                pl = [pid(x, yy) for yy in (y - 1, y) if 0 <= yy < Py]
                links.append(((x, y), (x + 1, y), pl))
        for x in range(Px + 1):
            for y in range(Py):                      # verticais
                pl = [pid(xx, y) for xx in (x - 1, x) if 0 <= xx < Px]
                links.append(((x, y), (x, y + 1), pl))
        self.links = links
        self.L = len(links)

    def link_bits(self, cfg: np.ndarray) -> np.ndarray:
        """cfg: inteiros de P bits (bit p = 1 ↔ μᶻ_p = −1). Retorna (n_cfg, L) bits de σˣ (1 ↔ σˣ = −1)."""
        mu = (cfg[:, None] >> np.arange(self.P)) & 1
        out = np.zeros((len(cfg), self.L), np.int64)
        for k, (_, _, pl) in enumerate(self.links):
            out[:, k] = mu[:, pl].sum(1) % 2
        return out

    def hamiltonian(self, h: float):
        n = 1 << self.P
        cfg = np.arange(n, dtype=np.int64)
        lb = self.link_bits(cfg)
        diag = -h * (1 - 2 * lb).sum(1).astype(float)
        rows, cols, vals = [cfg], [cfg], [diag]
        for p in range(self.P):
            rows.append(cfg); cols.append(cfg ^ (1 << p)); vals.append(np.full(n, -1.0))
        return coo_matrix((np.concatenate(vals), (np.concatenate(rows), np.concatenate(cols))), shape=(n, n)).tocsr()

    def ground_state(self, h: float):
        e, V = eigsh(self.hamiltonian(h), k=2, which="SA", tol=1e-12)
        o = np.argsort(e)
        psi = V[:, o[0]]
        return float(e[o[0]]), float(e[o[1]] - e[o[0]]), psi / np.linalg.norm(psi)

    def rect_links(self, x0, x1, y0, y1):
        """Elos que bordejam algum plaquete do retângulo [x0,x1)×[y0,y1)."""
        P = {x * self.Py + y for x in range(x0, x1) for y in range(y0, y1)}
        return [k for k, (_, _, pl) in enumerate(self.links) if P & set(pl)]

    def cut_vertices(self, A):
        Aset = set(A)
        inc = {}
        for k, (v1, v2, _) in enumerate(self.links):
            for v in (v1, v2):
                inc.setdefault(v, []).append(k)
        return [v for v, ks in inc.items() if any(k in Aset for k in ks) and any(k not in Aset for k in ks)], inc


def _H(p):
    p = np.asarray(p, float)
    p = p[p > 1e-15]
    return float(-(p * np.log(p)).sum())


def region_entropies(lat: Z2Lattice, psi: np.ndarray, A, lb: np.ndarray | None = None):
    """S_ext (Hilbert estendido), H_borda (centro elétrico) e S_dest = S_ext − H_borda."""
    if lb is None:
        lb = lat.link_bits(np.arange(len(psi), dtype=np.int64))
    Ac = [k for k in range(lat.L) if k not in set(A)]
    w = lambda cols: (lb[:, cols] << np.arange(len(cols))).sum(1)
    a, b = w(A), w(Ac)
    ua, ia = np.unique(a, return_inverse=True)
    ub, ib = np.unique(b, return_inverse=True)
    M = coo_matrix((psi, (ia, ib)), shape=(len(ua), len(ub))).tocsr()
    G = (M @ M.T) if len(ua) <= len(ub) else (M.T @ M)
    ev = np.linalg.eigvalsh(G.toarray())
    S_ext = _H(ev)
    cuts, inc = lat.cut_vertices(A)
    Aset = set(A)
    f = np.zeros(len(psi), np.int64)
    for j, v in enumerate(cuts):
        ks = [k for k in inc[v] if k in Aset]
        f |= (lb[:, ks].sum(1) % 2) << j
    pf = np.bincount(f, weights=psi ** 2)
    H_edge = _H(pf)
    return {"S_ext": S_ext, "H_edge": H_edge, "S_dist": S_ext - H_edge, "n_cut": len(cuts), "nA": len(A)}
