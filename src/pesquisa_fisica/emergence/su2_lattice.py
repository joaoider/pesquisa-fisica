"""Teoria de gauge SU(2) na rede em 1+1D com férmions (análogo não abeliano do Schwinger), exato.

Formulação hamiltoniana com férmions staggered de 2 cores (Bañuls, Cichy, Cirac, Jansen, Kühn,
PRX 7, 041046, 2017). Com condições abertas, os links são removidos por transformação de gauge e a
lei de Gauss não abeliana fixa o campo elétrico pela carga de cor acumulada:

    H = w Σ_{n,c} (φ†_{n,c} φ_{n+1,c} + h.c.) + m Σ_n (-1)^n n_n + J Σ_{n=0}^{N-2} S_n · S_n
    S_n = Σ_{k≤n} Q_k,   Q_k^a = φ†_k (σ^a/2) φ_k   (carga de cor no sítio k, a = x, y, z)

S_n · S_n é o Casimir do campo elétrico de cor na ligação n: acopla Q_k · Q_l para TODOS os k, l ≤ n,
inclusive termos de troca de cor φ†_{k0}φ_{k1}φ†_{l1}φ_{l0} entre sítios distantes (não abeliano).
Setor físico: N partículas (meio preenchimento) e singleto de cor total (S_tot² = 0), imposto por
penalidade P·S_tot².

Modos fermiônicos: q = 2n + c (bit q do índice da base = ocupação). Sinais de Jordan-Wigner explícitos.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.sparse.linalg import LinearOperator, eigsh


def _popcount(x):
    return np.bitwise_count(x).astype(np.int64)


class SU2Chain:
    def __init__(self, N: int, w=1.0, m=0.0, J=1.0, penalty=10.0):
        self.N, self.w, self.m, self.J, self.P = N, w, m, J, penalty
        self.n_modes = 2 * N
        self.basis = np.array(sorted(sum(1 << q for q in c) for c in combinations(range(2 * N), N)), dtype=np.int64)
        self.dim = len(self.basis)
        self._lookup = {int(s): i for i, s in enumerate(self.basis)}
        bits = (self.basis[:, None] >> np.arange(2 * N)) & 1
        n0, n1 = bits[:, 0::2], bits[:, 1::2]
        self.Qz = 0.5 * (n0 - n1)                                    # (dim, N)
        stag = (-1.0) ** np.arange(N)
        self.diag_mass = m * ((n0 + n1) * stag).sum(1)
        self.hops = [self._hop_map(2 * n + c, 2 * (n + 1) + c) for n in range(N - 1) for c in (0, 1)]
        self.qplus = [self._hop_map(2 * k, 2 * k + 1) for k in range(N)]   # Q⁺_k = φ†_{k0} φ_{k1}

    def _hop_map(self, a: int, b: int):
        """c†_a c_b (a < b): devolve (src, tgt, sinal) para os estados em que a ação é não nula."""
        s = self.basis
        ok = (((s >> b) & 1) == 1) & (((s >> a) & 1) == 0)
        src = np.nonzero(ok)[0]
        t = s[src] ^ ((1 << a) | (1 << b))
        between = ((1 << b) - 1) ^ ((1 << (a + 1)) - 1)
        sign = 1 - 2 * (_popcount(s[src] & between) & 1)
        tgt = np.array([self._lookup[int(x)] for x in t], dtype=np.int64)
        return src, tgt, sign.astype(float)

    @staticmethod
    def _apply(mp, v, dagger=False):
        src, tgt, sign = mp
        out = np.zeros_like(v)
        if not dagger:
            np.add.at(out, tgt, sign * v[src])
        else:
            np.add.at(out, src, sign * v[tgt])
        return out

    def casimir_apply(self, v, n: int):
        """(S_n · S_n) v = [S_z² + ½(S⁺S⁻ + S⁻S⁺)] v com S = Σ_{k≤n} Q_k."""
        Sz = self.Qz[:, : n + 1].sum(1)
        out = Sz**2 * v
        Sm_v = sum(self._apply(self.qplus[k], v, dagger=True) for k in range(n + 1))
        Sp_v = sum(self._apply(self.qplus[k], v) for k in range(n + 1))
        out += 0.5 * sum(self._apply(self.qplus[k], Sm_v) for k in range(n + 1))
        out += 0.5 * sum(self._apply(self.qplus[k], Sp_v, dagger=True) for k in range(n + 1))
        return out

    def matvec(self, v):
        v = np.asarray(v).ravel()
        out = self.diag_mass * v
        for mp in self.hops:
            out += self.w * (self._apply(mp, v) + self._apply(mp, v, dagger=True))
        for n in range(self.N - 1):
            out += self.J * self.casimir_apply(v, n)
        out += self.P * self.casimir_apply(v, self.N - 1)
        return out

    def ground_state(self, seed: int = 0):
        op = LinearOperator((self.dim, self.dim), matvec=self.matvec, dtype=float)
        v0 = np.random.default_rng(seed).standard_normal(self.dim)
        w, V = eigsh(op, k=1, which="SA", v0=v0, tol=1e-11)
        return float(w[0]), V[:, 0] / np.linalg.norm(V[:, 0])

    def total_color_casimir(self, v) -> float:
        return float(v @ self.casimir_apply(v, self.N - 1))

    def free_energy(self) -> float:
        """J = 0, sem penalidade: duas cópias (cores) de férmions livres staggered, N/2 partículas cada."""
        h = np.diag(self.m * (-1.0) ** np.arange(self.N))
        for n in range(self.N - 1):
            h[n, n + 1] = h[n + 1, n] = self.w
        e = np.sort(np.linalg.eigvalsh(h))
        return float(2 * e[: self.N // 2].sum())


# ---------------------------------------------------------------- emaranhamento com partição FERMIÔNICA

def fermionic_rdm(basis: np.ndarray, amps: np.ndarray, n_modes: int, modes: list[int]) -> np.ndarray:
    """Matriz densidade reduzida fermiônica dos modos `modes` (ordenados).

    Reordena os modos levando A para a frente (preservando a ordem interna); o sinal de cada estado é
    (-1)^{Σ_{a∈A} #(modos ocupados fora de A com índice < a)}. Depois do reordenamento, o traço parcial
    de qubits é o traço parcial fermiônico correto.
    """
    modes = sorted(modes)
    Aset = set(modes)
    rest = [q for q in range(n_modes) if q not in Aset]
    s = basis
    sign = np.zeros(len(s), dtype=np.int64)
    for a in modes:
        occ_a = (s >> a) & 1
        mask_before = sum(1 << q for q in rest if q < a)
        sign += occ_a * _popcount(s & mask_before)
    amp = amps * (1 - 2 * (sign & 1))
    x = sum(((s >> a) & 1) << i for i, a in enumerate(modes))
    r = sum(((s >> q) & 1) << i for i, q in enumerate(rest))
    M = np.zeros((1 << len(modes), 1 << len(rest)))
    np.add.at(M, (x, r), amp)
    return M @ M.T


def von_neumann(rho: np.ndarray) -> float:
    p = np.linalg.eigvalsh(rho)
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log(p)))


def site_mutual_information(basis, amps, n_modes: int, modes_per_site: int = 2) -> np.ndarray:
    n_sites = n_modes // modes_per_site
    site_modes = [list(range(modes_per_site * i, modes_per_site * (i + 1))) for i in range(n_sites)]
    s1 = [von_neumann(fermionic_rdm(basis, amps, n_modes, m)) for m in site_modes]
    I = np.zeros((n_sites, n_sites))
    for i in range(n_sites):
        for j in range(i + 1, n_sites):
            s12 = von_neumann(fermionic_rdm(basis, amps, n_modes, site_modes[i] + site_modes[j]))
            I[i, j] = I[j, i] = max(s1[i] + s1[j] - s12, 0.0)
    return I
