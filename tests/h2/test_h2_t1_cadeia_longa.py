"""H2-T1 cadeia longa (v0.14): razão assintótica das taxas de decaimento da IM, Schwinger m=0 vs férmion livre
a MESMO gap, via DMRG (TeNPy) até N=96.

Previsão pré-registrada (Notion, antes de rodar): razão → 2 ± 0,2 no contínuo (J→0).
P1: no gauge a correlação carregada decai MAIS RÁPIDO que a densidade (tensão de corda).
P2: no livre κ_D = 2 κ_C (Wick).
Resultado: P2 ✓; P1 ✗ (blindagem em m=0, não confinamento); razão → ≈1,5–1,6 ⇒ previsão de 2 FALSIFICADA.
"""
import json
from pathlib import Path

import numpy as np
import pytest

tenpy = pytest.importorskip("tenpy")
from scipy.sparse.linalg import eigsh  # noqa: E402

from pesquisa_fisica.emergence import schwinger as sw  # noqa: E402
from pesquisa_fisica.emergence.schwinger_dmrg import (  # noqa: E402
    decay_rate_power_fit, free_profiles, gauge_profiles, ground_and_first_excited)

RES = Path(__file__).resolve().parents[2] / "results"


def _rows():
    rows = []
    for fn in ["cadeia_longa_N40.json", "cadeia_longa_N64_96.json"]:
        rows += json.load(open(RES / fn))
    for r in rows:
        p = r["profiles"]
        d = np.array(p["d"])
        r["kg"] = decay_rate_power_fit(d, p["MI_f_gauge"])
        r["kf"] = 2 * np.arcsinh(r["m_free"] / 2)          # taxa assintótica exata da IM livre
        r["ratio"] = r["kg"] / r["kf"]
        r["kD_over_kC_g"] = decay_rate_power_fit(d, p["D_gauge"]) / decay_rate_power_fit(d, p["C_gauge"])
    return rows


@pytest.mark.parametrize("m,J", [(0.0, 0.0), (0.0, 1.0), (0.3, 2.0)])
def test_mpo_schwinger_reproduz_diagonalizacao_exata(m, J):
    N = 12
    H, _ = sw.hamiltonian(N, w=1.0, m=m, J=J)
    e = np.sort(eigsh(H, k=2, which="SA", tol=1e-12)[0])
    r = ground_and_first_excited(N, 1.0, m, J, chi_max=60)
    assert r["E0"] == pytest.approx(e[0], abs=1e-7)
    assert r["gap"] == pytest.approx(e[1] - e[0], abs=1e-5)


def test_referencia_livre_exata_igual_dmrg_em_J0():
    N, i, js = 20, 5, range(7, 15)
    g = gauge_profiles(N, 0.0, i, js, m=0.4, chi_max=60)
    f = free_profiles(N, 0.4, i, js)
    assert g["gap"] == pytest.approx(f["gap"], abs=1e-7)
    for k in ["C", "D", "MI_f"]:
        assert np.max(np.abs(g[k] - f[k])) < 1e-7


def test_P2_livre_densidade_decai_ao_dobro_da_carregada():
    for r in _rows():
        assert r["P2_kD_over_kC_free"] == pytest.approx(2.0, abs=0.07)  # N=96,J=1: D atinge o piso 1e-13


def test_P1_falhou_carregada_decai_mais_devagar_que_densidade_no_gauge():
    """Em m=0 o Schwinger BLINDA cargas: não há tensão de corda acelerando a correlação carregada."""
    for r in _rows():
        assert r["kD_over_kC_g"] > 1.3


def test_razao_convergida_em_N():
    rows = {(r["N"], r["J"]): r for r in _rows()}
    for J in [0.1, 0.25, 0.5, 1.0]:
        assert abs(rows[(64, J)]["ratio"] - rows[(96, J)]["ratio"]) < 0.06


def test_previsao_de_razao_2_falsificada():
    rows = [r for r in _rows() if r["N"] == 96]
    J = np.array([r["J"] for r in rows])
    R = np.array([r["ratio"] for r in rows])
    assert np.all((R > 1.35) & (R < 1.7))                   # gauge decai ~1,5× mais rápido que livre
    assert np.all(np.diff(R[np.argsort(J)]) < 0)            # cresce rumo ao acoplamento fraco
    small = J <= 0.25
    R0 = np.polyval(np.polyfit(J[small], R[small], 1), 0.0)  # extrapolação J → 0
    assert 1.45 < R0 < 1.75
    assert not (1.8 <= R0 <= 2.2)                           # critério de falha pré-registrado atingido
