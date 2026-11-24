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
