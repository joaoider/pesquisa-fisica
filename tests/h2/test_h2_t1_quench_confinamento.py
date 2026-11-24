"""H2-T1 confinamento na dinâmica (v0.16): quench global do vácuo nu, TDVP (N=24, t ≤ 6) vs livre exato de mesmo gap.

Pré-registro (Notion): P1 com m>0, S_g/S_f(t_max) < 0,5 e caindo com J; P2 em m=0 razão > 0,7;
P3 velocidade da frente de IM gauge < 0,6 da livre.
Resultado: P1 ✗ como enunciado (< 0,5 só em J=2; a tendência com J está certa); P2 ✓;
P3 ✗ (a frente gauge não tem velocidade: arranca e TRAVA, e trava também em m=0).
"""
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("tenpy")
from pesquisa_fisica.emergence.quench import free_quench, gauge_quench  # noqa: E402

RES = Path(__file__).resolve().parents[2] / "results"


def _rows(fn="quench_N24.json"):
    return {(r["m"], r["J"]): r for r in json.load(open(RES / fn))}


def _ratio(r):
    return r["S_gauge"][-1] / r["S_free"][-1]


def _reach(MI, thr=1e-3):
    MI = np.array(MI)
    hit = np.where((MI > thr).any(axis=0))[0]
    return int(hit.max()) + 1 if len(hit) else 0


def test_tdvp_reproduz_dinamica_livre_exata():
    times = np.arange(0, 2.01, 1.0)
    Sg, MIg, _ = gauge_quench(12, 0.5, 0.0, times, 3, [1, 2, 4], chi_max=32)
    Sf, MIf = free_quench(12, 0.5, times, 3, [1, 2, 4])
    assert np.max(np.abs(Sg - Sf)) < 1e-3 and np.max(np.abs(MIg - MIf)) < 1e-3


def test_P1_tendencia_certa_mas_limiar_falhou():
    R = _rows()
    r = [_ratio(R[(0.5, J)]) for J in (0.5, 1.0, 2.0)]
    assert r[0] > r[1] > r[2]                                  # supressão cresce com J (tensão de corda)
    assert r[2] < 0.5                                          # só J = 2 passa do limiar pré-registrado
    assert r[0] >= 0.5 and r[1] >= 0.5                         # → P1, como enunciado, FALHOU


def test_P2_blindagem_nao_suprime_e_massa_suprime_a_J_fixo():
    R = _rows()
    assert _ratio(R[(0.0, 0.5)]) > 0.7 and _ratio(R[(0.0, 1.0)]) > 0.7
    for J in (0.5, 1.0):                                       # mesma J: ligar m (corda) reduz o emaranhamento
        assert _ratio(R[(0.5, J)]) < _ratio(R[(0.0, J)]) - 0.1


def test_P3_falhou_frente_de_IM_trava_inclusive_sem_confinamento():
    R = _rows()
    for k, r in R.items():
        assert _reach(r["MI_gauge"]) < _reach(r["MI_free"])   # gauge: frente trava
    assert _reach(R[(0.0, 1.0)]["MI_gauge"]) == _reach(R[(0.5, 1.0)]["MI_gauge"])  # mas trava igual em m=0


def test_convergencia_em_chi():
    a, b = _rows(), _rows("quench_N24_chi160.json")
    for k in b:
        assert abs(_ratio(a[k]) - _ratio(b[k])) < 0.02
