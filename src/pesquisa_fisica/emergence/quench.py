"""Quench global no Schwinger (TDVP, TeNPy) e referência livre exata (gaussiana).

Estado inicial: vácuo nu staggered (sítios pares vazios, ímpares cheios — Q_k = 0, sem pares, sem fluxo).
Observáveis: entropia de meia cadeia S(t) e IM de dois sítios (partição fermiônica) I(c, c+d; t).
"""
from __future__ import annotations

import numpy as np
from scipy.linalg import eigh, expm
from tenpy.algorithms.tdvp import TwoSiteTDVPEngine
from tenpy.networks.mps import MPS

from .schwinger_dmrg import schwinger_mpo_model, two_site_mi


def _vn(p):
    p = np.clip(np.asarray(p, float), 1e-15, 1 - 1e-15)
    return float(-(p * np.log(p) + (1 - p) * np.log(1 - p)).sum())


def free_quench(N: int, m: float, times, center: int, ds, w: float = 1.0):
    h = np.diag([m * (-1.0) ** n for n in range(N)]) + w * (np.eye(N, k=1) + np.eye(N, k=-1))
    C0 = np.diag([0.0 if n % 2 == 0 else 1.0 for n in range(N)])
    e, U = eigh(h)
    S, MI = [], []
    for t in times:
        Ut = U @ np.diag(np.exp(-1j * e * t)) @ U.conj().T
        C = Ut.conj() @ C0 @ Ut.T                        # C_ij = ⟨c_i† c_j⟩(t)
        S.append(_vn(np.linalg.eigvalsh(C[: N // 2, : N // 2])))
        n = np.real(np.diag(C))
        row = []
        for d in ds:
            j = center + d
            x = abs(C[center, j])
            row.append(two_site_mi(n[center], n[j], n[center] * n[j] - x ** 2, x))
        MI.append(row)
    return np.array(S), np.array(MI)


def gauge_quench(N: int, m: float, J: float, times, center: int, ds, w: float = 1.0,
                 chi_max: int = 128, dt: float = 0.05):
    M = schwinger_mpo_model(N, w, m, J)
    init = ["down" if k % 2 == 0 else "up" for k in range(N)]
    psi = MPS.from_product_state(M.lat.mps_sites(), init, bc="finite")
    eng = TwoSiteTDVPEngine(psi, M, {"dt": dt, "N_steps": 1,
                                      "trunc_params": {"chi_max": chi_max, "svd_min": 1e-10}})
    js = [center + d for d in ds]
    S, MI, chis = [], [], []
    t_now = 0.0
    for t in times:
        while t_now < t - 1e-9:
            eng.run()
            t_now += dt
        S.append(float(psi.entanglement_entropy()[N // 2 - 1]))
        n = np.real(psi.expectation_value("Sz")) + 0.5
        zz = np.real(psi.correlation_function("Sz", "Sz", sites1=[center], sites2=js)[0])
        Cf = np.abs(psi.correlation_function("Sp", "Sm", sites1=[center], sites2=js, opstr="Sigmaz",
                                             str_on_first=False)[0])
        ninj = zz + 0.5 * (n[center] - 0.5) + 0.5 * (n[js] - 0.5) + 0.25
        MI.append([two_site_mi(n[center], n[j], float(ninj[k]), float(Cf[k])) for k, j in enumerate(js)])
        chis.append(max(psi.chi))
    return np.array(S), np.array(MI), chis


# ---------------------------------------------------------------- MDRs de blocos disjuntos (partição de spin)

def block_pair_rdms(psi: MPS, a0: int, b_starts, size: int = 2) -> dict:
    """ρ(A ∪ B) para A = [a0, a0+size) fixo e B = [b, b+size) para cada b em b_starts (b ≥ a0+size).
    Contração direta do MPS na forma B (canônica à direita); custo O(N χ³ 4^size)."""
    b_starts = sorted(b_starts)
    S = psi.get_SL(a0)
    env = np.diag(np.asarray(S) ** 2)[None, :, :].astype(complex)          # (K, χ, χ)
    B = lambda i: psi.get_B(i, form="B").to_ndarray()                       # (vL, p, vR)

    def open_site(env, i):
        T = B(i)
        tmp = np.tensordot(env, T, axes=([1], [0]))                      # k b p c
        new = np.tensordot(tmp, T.conj(), axes=([1], [0]))                # k p c q d
        new = new.transpose(0, 1, 3, 2, 4)                                # k p q c d
        K, d = env.shape[0], T.shape[1]
        return new.reshape(K * d * d, T.shape[2], T.shape[2])

    def close_site(env, i):
        T = B(i)
        tmp = np.tensordot(env, T, axes=([1], [0]))                      # k b p c
        return np.tensordot(tmp, T.conj(), axes=([1, 2], [0, 1]))        # k c d

    for i in range(a0, a0 + size):
        env = open_site(env, i)
    out = {}
    pos = a0 + size
    for b in b_starts:
        while pos < b:
            env = close_site(env, pos)
            pos += 1
        e = env
        for i in range(b, b + size):
            e = open_site(e, i)
        vec = np.einsum("kcc->k", e)                                         # traço à direita
        n = 2 * size                                                          # sítios abertos
        t = vec.reshape([2, 2] * n)                                          # (p0,q0,p1,q1,...)
        t = t.transpose(list(range(0, 2 * n, 2)) + list(range(1, 2 * n, 2)))
        out[b] = t.reshape(2 ** n, 2 ** n)
    return out


def _S_rho(r):
    w = np.linalg.eigvalsh((r + r.conj().T) / 2)
    w = w[w > 1e-14]
    return float(-(w * np.log(w)).sum())


def _ptrace(rho, keep, n):
    t = rho.reshape([2] * (2 * n))
    drop = [i for i in range(n) if i not in keep]
    for k, i in enumerate(sorted(drop, reverse=True)):
        m = n - k
        t = np.trace(t, axis1=i, axis2=i + m)
    d = 2 ** len(keep)
    return t.reshape(d, d)


def block_mi_and_excess(rho4):
    """Para ρ sobre (a1,a2,b1,b2): I(A:B) e Σ I(a:b) dos 4 pares de sítios."""
    n = 4
    SA, SB, SAB = _S_rho(_ptrace(rho4, [0, 1], n)), _S_rho(_ptrace(rho4, [2, 3], n)), _S_rho(rho4)
    S1 = [_S_rho(_ptrace(rho4, [i], n)) for i in range(n)]
    pair = 0.0
    for a in (0, 1):
        for b in (2, 3):
            pair += S1[a] + S1[b] - _S_rho(_ptrace(rho4, [a, b], n))
    return SA + SB - SAB, pair


def block_entropy(psi: MPS, a: int, ell: int) -> float:
    """S do bloco contíguo [a, a+ell) via SVD de X[(p…),(α,γ)] = S_α (B…B)_{α p… γ}. Dimensão 2^ℓ × χ_L χ_R."""
    S = np.asarray(psi.get_SL(a))
    T = np.diag(S).astype(complex)                                   # (α, χ)
    T = T[:, None, :]                                                 # (α, P=1, χ)
    for i in range(a, a + ell):
        Bi = psi.get_B(i, form="B").to_ndarray()                      # (χ, p, χ')
        T = np.tensordot(T, Bi, axes=([2], [0]))                      # (α, P, p, χ')
        T = T.reshape(T.shape[0], T.shape[1] * T.shape[2], T.shape[3])
    X = T.transpose(1, 0, 2).reshape(T.shape[1], -1)
    s = np.linalg.svd(X, compute_uv=False) ** 2
    s = s[s > 1e-14]
    return float(-(s * np.log(s)).sum())


def cut_entropy(psi: MPS, x: int) -> float:
    """Entropia do corte que separa [0, x) de [x, N)."""
    return float(psi.entanglement_entropy()[x - 1])


def tripartite_info(psi: MPS, a: int, size: int = 2) -> float:
    """I₃(A:B:C) para blocos adjacentes A=[a,a+s), B=[a+s,a+2s), C=[a+2s,a+3s)."""
    s = size
    SA, SB, SC = (block_entropy(psi, a + k * s, s) for k in range(3))
    SAB, SBC = block_entropy(psi, a, 2 * s), block_entropy(psi, a + s, 2 * s)
    SABC = block_entropy(psi, a, 3 * s)
    SAC = _S_rho(block_pair_rdms(psi, a, [a + 2 * s], size=s)[a + 2 * s])
    return SA + SB + SC - SAB - SBC - SAC + SABC
