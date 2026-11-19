"""Teoria de gauge SU(N_c) na rede em 1+1D com férmions staggered — geral (SU(2), SU(3) = QCD 1+1D).

Mesma construção de su2_lattice.py (links removidos por transformação de gauge, lei de Gauss não abeliana
fixa o campo elétrico de cor pela carga acumulada), para N_c cores:

    H = w Σ_{n,c} (φ†_{n,c} φ_{n+1,c} + h.c.) + m Σ_n (-1)^n n_n + J Σ_{n=0}^{N-2} S_n·S_n + P·S_tot·S_tot

Casimir pela identidade de Fierz  Σ_a T^a_{ij} T^a_{kl} = ½(δ_il δ_jk − δ_ij δ_kl / N_c):
    S_n·S_n = ½ Σ_{c,c'} E^{cc'}_n E^{c'c}_n − (1 / 2N_c) (N_n)²,
    E^{cc'}_n = Σ_{k≤n} φ†_{k,c} φ_{k,c'}   (troca de cor; E^{cc}_n = nº de partículas de cor c até n),
    N_n = Σ_{k≤n} n_k.
Setor: meio preenchimento (N·N_c/2 partículas) e singleto de cor total (penalidade P).
Modos: q = N_c·n + c. Mapas de salto são injetivos (sem colisões), então a aplicação é vetorizada.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


def _popcount(x):
    return np.bitwise_count(x).astype(np.int64)


def half_filling_basis(n_modes: int, n_particles: int) -> np.ndarray:
    if n_modes <= 26:
        allst = np.arange(1 << n_modes, dtype=np.int64)
        return allst[_popcount(allst) == n_particles]
    from itertools import combinations
    return np.array(sorted(sum(1 << q for q in c) for c in combinations(range(n_modes), n_particles)), dtype=np.int64)


class SUNChain:
    def __init__(self, N: int, Nc: int = 3, w=1.0, m=0.0, J=1.0, penalty=10.0):
        if (N * Nc) % 2:
            raise ValueError("N·N_c precisa ser par (meio preenchimento)")
        self.N, self.Nc, self.w, self.m, self.J, self.P = N, Nc, w, m, J, penalty
        self.n_modes = N * Nc
        self.basis = half_filling_basis(self.n_modes, N * Nc // 2)
        self.dim = len(self.basis)
        bits = ((self.basis[:, None] >> np.arange(self.n_modes)) & 1).astype(np.int8)
        occ = bits.reshape(self.dim, N, Nc)                                  # occ[:, site, color]
        self.cum_color = np.cumsum(occ, axis=1).astype(np.int8)              # (dim, N, Nc): nº de cor c até o sítio n
        stag = (-1.0) ** np.arange(N)
        self.diag_mass = m * (occ.sum(2) * stag).sum(1)
        del bits, occ
        self.hops = [self._map(Nc * n + c, Nc * (n + 1) + c) for n in range(N - 1) for c in range(Nc)]
        # E^{cc'}_k para c ≠ c' (troca de cor no sítio k)
        self.E = {(c, cp): [self._map(Nc * k + c, Nc * k + cp) for k in range(N)]
                  for c in range(Nc) for cp in range(Nc) if c != cp}

    def _map(self, a: int, b: int):
        """c†_a c_b (a ≠ b): (src, tgt, sinal) nos estados em que a ação é não nula. Injetivo."""
        s = self.basis
        ok = (((s >> b) & 1) == 1) & (((s >> a) & 1) == 0)
        src = np.nonzero(ok)[0]
        t = s[src] ^ ((1 << a) | (1 << b))
        lo, hi = min(a, b), max(a, b)
        between = ((1 << hi) - 1) ^ ((1 << (lo + 1)) - 1)
        sign = (1 - 2 * (_popcount(s[src] & between) & 1)).astype(np.int8)
        tgt = np.searchsorted(s, t)
        return src.astype(np.int32), tgt.astype(np.int32), sign

    @staticmethod
    def _apply(mp, v, out=None, scale=1.0):
        src, tgt, sign = mp
        if out is None:
            out = np.zeros_like(v)
        out[tgt] += scale * sign * v[src]
        return out

    def casimir_apply(self, v, n: int):
        """(S_n·S_n) v."""
        cum = self.cum_color[:, n, :].astype(np.float64)
        Ntot = cum.sum(1)
        out = (0.5 * (cum**2).sum(1) - Ntot**2 / (2 * self.Nc)) * v
        for (c, cp), maps in self.E.items():
            inner = np.zeros_like(v)
            for k in range(n + 1):
                self._apply(maps[k], v, inner)              # E^{c'c}... usamos o par (c, cp) e depois (cp, c)
            outer_maps = self.E[(cp, c)]
            for k in range(n + 1):
                self._apply(outer_maps[k], inner, out, 0.5)
        return out

    def _link_coefficients(self):
        c = np.full(self.N, self.J, dtype=float)
        c[-1] = self.P                     # última "ligação" = carga de cor total (penalidade de singleto)
        return c

    def matvec(self, v):
        """Hv com o Casimir de todas as ligações calculado de forma incremental:
        E_{≤n} v = E_{≤n-1} v + E_n v (uma aplicação por n para o fator interno)."""
        v = np.asarray(v, dtype=np.float64).ravel()
        out = self.diag_mass * v
        for src, tgt, sign in self.hops:
            out[tgt] += self.w * sign * v[src]
            out[src] += self.w * sign * v[tgt]
        coef = self._link_coefficients()
        # parte diagonal de todas as ligações
        diag = np.zeros(self.dim)
        for n in range(self.N):
            cum = self.cum_color[:, n, :].astype(np.float64)
            diag += coef[n] * (0.5 * (cum**2).sum(1) - cum.sum(1) ** 2 / (2 * self.Nc))
        out += diag * v
        # parte de troca de cor: ½ Σ_{c≠c'} E^{c'c}_{≤n} E^{cc'}_{≤n}
        for (c, cp), maps in self.E.items():
            outer_maps = self.E[(cp, c)]
            inner = np.zeros(self.dim)
            for n in range(self.N):
                self._apply(maps[n], v, inner)
                for k in range(n + 1):
                    self._apply(outer_maps[k], inner, out, 0.5 * coef[n])
        return out

    def ground_state(self, seed: int = 0, tol: float = 1e-10):
        op = LinearOperator((self.dim, self.dim), matvec=self.matvec, dtype=np.float64)
        v0 = np.random.default_rng(seed).standard_normal(self.dim)
        w, V = eigsh(op, k=1, which="SA", v0=v0, tol=tol)
        return float(w[0]), V[:, 0] / np.linalg.norm(V[:, 0])

    def total_color_casimir(self, v) -> float:
        return float(v @ self.casimir_apply(v, self.N - 1))

    def free_energy(self) -> float:
        """J = 0: N_c cópias de férmions livres staggered, N/2 partículas cada."""
        h = np.diag(self.m * (-1.0) ** np.arange(self.N))
        for n in range(self.N - 1):
            h[n, n + 1] = h[n + 1, n] = self.w
        e = np.sort(np.linalg.eigvalsh(h))
        return float(self.Nc * e[: self.N // 2].sum())


# ---------------------------------------------------------------- emaranhamento (partição fermiônica)

def fermionic_rdm(basis, amps, n_modes: int, modes) -> np.ndarray:
    modes = sorted(modes)
    Aset = set(modes)
    rest = [q for q in range(n_modes) if q not in Aset]
    s = basis
    parity = np.zeros(len(s), dtype=np.int64)
    for a in modes:
        mask_before = sum(1 << q for q in rest if q < a)
        parity += ((s >> a) & 1) * _popcount(s & mask_before)
    amp = amps * (1 - 2 * (parity & 1))
    x = np.zeros(len(s), dtype=np.int64)
    for i, a in enumerate(modes):
        x |= ((s >> a) & 1) << i
    r = np.zeros(len(s), dtype=np.int64)
    for i, q in enumerate(rest):
        r |= ((s >> q) & 1) << i
    # compacta as colunas usadas (evita matriz 2^|rest| inteira)
    r_unique, r_idx = np.unique(r, return_inverse=True)
    M = np.zeros((1 << len(modes), len(r_unique)))
    M[x, r_idx] = amp
    return M @ M.T


def von_neumann(rho) -> float:
    p = np.linalg.eigvalsh(rho)
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log(p)))


def site_mutual_information(chain: SUNChain, amps) -> np.ndarray:
    Nc, N = chain.Nc, chain.N
    sm = [list(range(Nc * i, Nc * (i + 1))) for i in range(N)]
    s1 = [von_neumann(fermionic_rdm(chain.basis, amps, chain.n_modes, m)) for m in sm]
    I = np.zeros((N, N))
    for i in range(N):
        for j in range(i + 1, N):
            s12 = von_neumann(fermionic_rdm(chain.basis, amps, chain.n_modes, sm[i] + sm[j]))
            I[i, j] = I[j, i] = max(s1[i] + s1[j] - s12, 0.0)
    return I
