"""H2-T1 coarse-graining em escala de méson (v0.19): quench EXATO (N=20, setor Q=0), IM entre blocos simétricos.

Pré-registro: P1 (m=0,J=1,ℓ_b=4) pico balístico em todo r e t_pico crescente; P2 v_gauge/v_livre ∈ [0,5;1];
P3 (m=0,5,J=2) I_max(r=8) < 0,1·livre.  Resultado: P1 ✗ (picos simultâneos t≈1,6 p/ todo r), P2 ✗, P3 ✓.
Diagnóstico: a IM total (partição de qubits, Gauss eliminada) é dominada por coerências de CORDA, com frente
sublinear e "rápida" (v≈3,8 > 2w). Só com cargas (observáveis locais invariantes de gauge) o cone de luz
é causal: em m=0 coincide com o livre (v≈2); com confinamento é freado.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from pesquisa_fisica.emergence.ed_quench import (_S, arrival_time, block_mi, classical_block_mi, evolve, rdm,
                                                 symmetric_blocks)
from pesquisa_fisica.emergence.quench import free_quench

RES = Path(__file__).resolve().parents[2] / "results"
RS = (2, 4, 6, 8)


def _main():
    return {(r["m"], r["J"]): r for r in json.load(open(RES / "meson_coarse_N20.json"))}


def _class():
    d = json.load(open(RES / "meson_coarse_classico_N20.json"))
    return d["times"], {(r["m"], r["J"], r["side"]): r for r in d["rows"]}


def test_evolucao_exata_reproduz_gaussiana():
    N, ts = 10, [0.0, 1.0, 2.0]
    st = evolve(N, 0.5, 0.0, ts)
    Sf, _ = free_quench(N, 0.5, ts, 2, [1])
    assert np.allclose([_S(rdm(s, N, list(range(N // 2)))) for s in st], Sf, atol=1e-9)


def test_blocos_simetricos():
    A, B = symmetric_blocks(20, 4, 8)
    assert A == [2, 3, 4, 5] and B == [14, 15, 16, 17]


def test_P1_P2_falharam_picos_simultaneos_no_gauge():
    R = _main()[(0.0, 1.0)]
    g = [R["gauge"][f"4_{r}"] for r in RS]
    assert not all(x["peak"] for x in g)                                  # sem pico balístico em r=6,8
    tp = np.array([x["t_peak"] for x in g])
    assert np.ptp(tp) < 0.1                                              # t_pico ≈ 1,6 para todo r
    tf = np.array([R["free"][f"4_{r}"]["t_peak"] for r in RS])
    assert np.all(np.diff(tf) > 0)                                       # livre: balístico


def test_P3_passou_confinamento_mata_correlacao_longa():
    R = _main()[(0.5, 2.0)]
    assert R["gauge"]["4_8"]["I_max"] / R["free"]["4_8"]["I_max"] < 0.1


def test_frente_total_do_gauge_e_rapida_e_sublinear():
    R = _main()[(0.0, 1.0)]
    t = R["times"]
    ta = np.array([arrival_time(t, R["gauge"][f"4_{r}"]["I"], 1e-3) for r in RS])
    v = 1 / (2 * np.polyfit(RS, ta, 1)[0])
    assert v > 3.0                                                       # acima de 2w = 2
    assert np.all(np.diff(np.diff(ta)) < 0)                              # incrementos decrescentes


def test_cone_de_luz_causal_nas_cargas():
    t, C = _class()
    arr = lambda k: np.array([arrival_time(t, C[k]["r"][str(r)]["I_class"], 1e-4) for r in RS])
    g0, f0 = arr((0.0, 1.0, "gauge")), arr((0.0, 1.0, "free"))
    assert np.max(np.abs(g0 - f0)) < 0.15                                # m=0: mesmo cone de luz do livre
    v0 = 1 / (2 * np.polyfit(RS, g0, 1)[0])
    assert 1.5 < v0 < 2.2
    g2, f2 = arr((0.5, 2.0, "gauge")), arr((0.5, 2.0, "free"))
    assert g2[-1] > f2[-1] + 1.0                                         # confinamento freia (r=8)


def test_parte_coerente_domina_a_IM_total_do_gauge():
    t, C = _class()
    rec = C[(0.0, 1.0, "gauge")]["r"]["8"]
    k = int(np.argmax(rec["I_full"]))
    assert rec["I_class"][k] < 0.05 * rec["I_full"][k]
