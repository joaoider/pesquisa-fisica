"""H2-T1 regime confinante (v0.15): Schwinger m ≠ 0 vs livre a mesmo gap (DMRG, N=64; dados em results/).

Pré-registro (Notion): P1 taxa local κ_C(d) da correlação carregada cresce com d (corda); P2 satura em d* ≈ gap/J;
P3 razão de IM gauge/livre sobe acima de 1,5 com m.  Resultado: P1 ✗, P2 ✗ (nada para saturar), P3 ✗ (CAI).
Achado: a razão é função só de x = m/√J (≡ m/g no contínuo): ≈1,5 em x=0 → 1 quando x → ∞.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from pesquisa_fisica.emergence.schwinger_dmrg import decay_rate_power_fit

RES = Path(__file__).resolve().parents[2] / "results"


def _rows():
    rows = json.load(open(RES / "confinante_N64.json")) + json.load(open(RES / "confinante_pesado_N64.json"))
    for r in rows:
        d = np.array(r["d"])
        r["ratio"] = decay_rate_power_fit(d, r["MI_f_gauge"]) / (2 * np.arcsinh(r["m_free"] / 2))
        r["x"] = r["m"] / np.sqrt(r["J"])
        C = np.array(r["C_gauge"])
        loc = {}
        for par in (0, 1):
            sel = (d % 2 == par) & (C > 1e-12)
            loc[par] = -np.diff(np.log(C[sel])) / 2
        r["loc"] = loc
    return rows


def test_P1_falhou_sem_assinatura_de_corda_no_correlador_estatico():
    """Critério pré-registrado: P1 falha se κ_C(d) for constante (10%) em toda a janela para m ≥ 0,5.
    Onde varia mais que 10% (quarks pesados), varia no sentido OPOSTO (κ cai); a maior subida é < 3%. Atingido em todos os pontos; com quarks pesados (corda de d* ≈ 21 sítios) κ_C(d) chega a DECRESCER."""
    rows = _rows()
    for r in rows:
        if r["m"] >= 0.5:
            for loc in r["loc"].values():
                rise = max(loc[k] - loc[:k + 1].min() for k in range(len(loc)))   # maior SUBIDA de κ_C(d)
                assert rise / loc.mean() < 0.03, (r["m"], r["J"])               # corda exigiria subida forte
    heavy = [r for r in rows if r["J"] == 0.1][0]
    for loc in heavy["loc"].values():
        assert loc[-1] < loc[0]


def test_P3_falhou_razao_cai_com_a_massa():
    rows = sorted([r for r in _rows() if r["J"] == 0.5], key=lambda r: r["m"])
    R = [r["ratio"] for r in rows]
    assert R[0] == pytest.approx(1.49, abs=0.05)
    assert np.all(np.diff(R) < 0) and R[-1] < 1.25


def test_razao_e_funcao_universal_de_m_sobre_g():
    rows = sorted(_rows(), key=lambda r: r["x"])
    x = np.array([r["x"] for r in rows]); R = np.array([r["ratio"] for r in rows])
    # pontos com J diferentes e mesmo x colapsam
    a = [r for r in rows if r["m"] == 0.25 and r["J"] == 0.5][0]
    b = [r for r in rows if r["m"] == 0.5 and r["J"] == 2.0][0]
    assert abs(a["x"] - b["x"]) < 1e-9 and abs(a["ratio"] - b["ratio"]) < 0.03
    # monótona em x (tolerância para pontos de mesmo x)
    order = np.argsort(x)
    assert np.all(np.diff(R[order]) < 0.02)
    assert R[x > 3].max() < 1.1                                             # quarks pesados → livre
