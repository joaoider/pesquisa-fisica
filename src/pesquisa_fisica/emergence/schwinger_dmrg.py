"""Modelo de Schwinger em cadeias longas via DMRG (TeNPy), mesma Hamiltoniana de schwinger.py.

H = w Σ (σ⁺_n σ⁻_{n+1} + h.c.) + (m/2) Σ (−1)^n σᶻ_n + J Σ_{n=0}^{N−2} L_n²,   L_n = Σ_{k≤n} Q_k,
Q_k = (σᶻ_k + (−1)^k)/2 = S^z_k + (−1)^k/2.

J Σ_n L_n² = J Σ_k (N−1−k) Q_k²  +  2J Σ_{k<l} (N−1−l) Q_k Q_l   (termos com índice ≤ N−2).
MPO de dimensão 5: canais [0 início, 1 σ⁺ pendente, 2 σ⁻ pendente, 3 soma acumulada de Q, 4 fim].
Q_k² = ½ + (−1)^k S^z_k (porque (S^z)² = ¼).
"""
from __future__ import annotations

import numpy as np
from tenpy.algorithms import dmrg
from tenpy.models.lattice import Chain
from tenpy.models.model import MPOModel
from tenpy.networks.mpo import MPO
from tenpy.networks.mps import MPS
from tenpy.networks.site import SpinHalfSite


def schwinger_mpo_model(N: int, w: float = 1.0, m: float = 0.0, J: float = 1.0) -> MPOModel:
    site = SpinHalfSite(conserve="Sz")
    lat = Chain(N, site, bc="open", bc_MPS="finite")
    grids = []
    for l in range(N):
        s = (-1.0) ** l
        cJ = J * (N - 1 - l) if l <= N - 2 else 0.0
        onsite = [("Sz", m * s), ("Sz", cJ * s), ("Id", 0.5 * cJ)]   # (m/2)σᶻ = m S^z ; J(N−1−l)Q²
        Q = [("Sz", 1.0), ("Id", 0.5 * s)]
        Qc = [("Sz", 2 * cJ), ("Id", cJ * s)]                        # 2J(N−1−l) Q_l
        g = [[None] * 5 for _ in range(5)]
        g[0][0] = "Id"
        g[4][4] = "Id"
        g[0][1] = [("Sp", w)]
        g[1][4] = "Sm"
        g[0][2] = [("Sm", w)]
        g[2][4] = "Sp"
        g[0][3] = Q
        g[3][3] = "Id"
        g[3][4] = Qc
        g[0][4] = onsite
        if l == 0:
            g = [g[0]]
        if l == N - 1:
            g = [[row[4]] for row in g]
        grids.append(g)
    H = MPO.from_grids(lat.mps_sites(), grids, bc="finite", IdL=0, IdR=-1)
    return MPOModel(lat, H)


def ground_and_first_excited(N: int, w=1.0, m=0.0, J=1.0, chi_max: int = 120, excited: bool = True):
    M = schwinger_mpo_model(N, w, m, J)
    init = ["down" if k % 2 == 0 else "up" for k in range(N)]          # vácuo staggered: Q_k = 0, S^z_tot = 0
    params = {"mixer": True, "max_E_err": 1e-12, "max_sweeps": 40,
              "trunc_params": {"chi_max": chi_max, "svd_min": 1e-14}}
    psi0 = MPS.from_product_state(M.lat.mps_sites(), init, bc="finite")
    E0 = dmrg.run(psi0, M, params)["E"]
    out = {"E0": E0, "psi0": psi0}
    if excited:
        init1 = list(init)
        init1[N // 2], init1[N // 2 + 1] = init1[N // 2 + 1], init1[N // 2]   # mesmo setor S^z = 0
        psi1 = MPS.from_product_state(M.lat.mps_sites(), init1, bc="finite")
        E1 = dmrg.run(psi1, M, params, orthogonal_to=[psi0.copy()])["E"]
        out.update(E1=E1, gap=E1 - E0)
    return out


def charged_correlator(psi: MPS, i: int, js) -> np.ndarray:
    """|⟨c_i† c_j⟩| = |⟨σ⁺_i (Π_{i<k<j} (−σᶻ_k)) σ⁻_j⟩| (Jordan-Wigner) — correlação CARREGADA."""
    vals = psi.correlation_function("Sp", "Sm", sites1=[i], sites2=list(js), opstr="Sigmaz", str_on_first=False)
    return np.abs(vals[0])


def density_correlator(psi: MPS, i: int, js) -> np.ndarray:
    """|⟨Q_i Q_j⟩_c| = |⟨S^z_i S^z_j⟩ − ⟨S^z_i⟩⟨S^z_j⟩| — correlação NEUTRA."""
    zz = psi.correlation_function("Sz", "Sz", sites1=[i], sites2=list(js))[0]
    sz = psi.expectation_value("Sz")
    return np.abs(zz - sz[i] * sz[list(js)])


# ---------------------------------------------------------------- cadeia livre exata (gaussiana)

def free_chain(N: int, w: float = 1.0, m: float = 0.0):
    """Férmion livre staggered (J = 0) exato: h = w(c†_n c_{n+1} + h.c.) + m(−1)^n n_n, meia ocupação.
    Devolve (C, gap) com C_ij = ⟨c_i† c_j⟩ e gap = ε_{N/2} − ε_{N/2−1} (partícula-buraco, mesmo setor de carga)."""
    h = np.diag([m * (-1.0) ** n for n in range(N)]) + w * (np.eye(N, k=1) + np.eye(N, k=-1))
    e, U = np.linalg.eigh(h)
    occ = U[:, : N // 2]
    return occ @ occ.T, e[N // 2] - e[N // 2 - 1]


def free_gap(N: int, m: float, w: float = 1.0) -> float:
    return free_chain(N, w, m)[1]


# ---------------------------------------------------------------- informação mútua de 2 sítios

def _H(p):
    p = np.asarray(p, float)
    p = p[p > 1e-15]
    return float(-(p * np.log(p)).sum())


def two_site_mi(ni, nj, ninj, off):
    """IM entre os modos i e j a partir de ⟨n_i⟩, ⟨n_j⟩, ⟨n_i n_j⟩ e do elemento fora da diagonal |x|.
    x = ⟨c_i† c_j⟩ (com string JW) → partição FERMIÔNICA;  x = ⟨σ⁺_i σ⁻_j⟩ (sem string) → partição de SPIN."""
    a, b, c = ninj, ni - ninj, nj - ninj
    d = 1 - ni - nj + ninj
    tr = (b + c) / 2
    disc = np.sqrt(((b - c) / 2) ** 2 + off ** 2)
    S_ij = _H([a, d, tr + disc, tr - disc])
    return _H([ni, 1 - ni]) + _H([nj, 1 - nj]) - S_ij


def free_profiles(N: int, m: float, i: int, js, w: float = 1.0):
    C, gap = free_chain(N, w, m)
    js = list(js)
    n = np.diag(C)
    Cij = C[i, js]
    D = Cij ** 2                                     # |⟨n_i n_j⟩_c| = |C_ij|² (Wick)
    mi_f = np.array([two_site_mi(n[i], n[j], n[i] * n[j] - C[i, j] ** 2, abs(C[i, j])) for j in js])
    return {"gap": gap, "C": np.abs(Cij), "D": D, "MI_f": mi_f}


def gauge_profiles(N: int, J: float, i: int, js, m: float = 0.0, w: float = 1.0, chi_max: int = 120):
    r = ground_and_first_excited(N, w, m, J, chi_max=chi_max)
    psi = r["psi0"]
    js = list(js)
    n = psi.expectation_value("Sz") + 0.5
    zz = psi.correlation_function("Sz", "Sz", sites1=[i], sites2=js)[0]
    ninj = zz + 0.5 * (n[i] - 0.5) + 0.5 * (n[js] - 0.5) + 0.25        # ⟨n_i n_j⟩ com n = S^z + ½
    Cf = charged_correlator(psi, i, js)
    Cs = np.abs(psi.correlation_function("Sp", "Sm", sites1=[i], sites2=js)[0])
    mi_f = np.array([two_site_mi(n[i], n[j], ninj[k], Cf[k]) for k, j in enumerate(js)])
    mi_s = np.array([two_site_mi(n[i], n[j], ninj[k], Cs[k]) for k, j in enumerate(js)])
    return {"gap": r["gap"], "E0": r["E0"], "C": Cf, "D": np.abs(ninj - n[i] * n[js]),
            "MI_f": mi_f, "MI_s": mi_s}


def match_free_mass(N: int, target_gap: float, w: float = 1.0, lo: float = 0.0, hi: float = 6.0, iters: int = 60):
    for _ in range(iters):
        mid = 0.5 * (lo + hi)
        if free_gap(N, mid, w) < target_gap:
            lo = mid
        else:
            hi = mid
    return 0.5 * (lo + hi)


def decay_rate_fit(ds, prof, floor: float = 1e-13, parity: int | None = 0):
    """Taxa κ de |prof| ~ e^{−κ d} por mínimos quadrados em log, só pontos acima do piso (e da paridade dada)."""
    ds = np.asarray(ds)
    prof = np.asarray(prof)
    sel = prof > floor
    if parity is not None:
        sel &= (ds % 2 == parity)
    if sel.sum() < 3:
        return np.nan
    return float(-np.polyfit(ds[sel], np.log(prof[sel]), 1)[0])


def decay_rate_power_fit(ds, prof, floor: float = 1e-13, parity: int | None = 0, d_min: int = 4):
    """κ de |prof| ~ A d^{−α} e^{−κ d} (ajuste linear em log com 3 parâmetros) — remove o viés do prefator."""
    ds = np.asarray(ds, float)
    prof = np.asarray(prof)
    sel = (prof > floor) & (ds >= d_min)
    if parity is not None:
        sel &= (ds.astype(int) % 2 == parity)
    if sel.sum() < 4:
        return np.nan
    X = np.column_stack([np.ones(sel.sum()), -ds[sel], -np.log(ds[sel])])
    coef, *_ = np.linalg.lstsq(X, np.log(prof[sel]), rcond=None)
    return float(coef[1])
