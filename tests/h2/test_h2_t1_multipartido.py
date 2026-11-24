"""H2-T1 correlações multipartidas após o quench (v0.17): N=24, t=6, gauge vs livre (ambos TDVP, partição de spin).

Pré-registro: P1 ganho de alcance blocos−sítios do gauge > livre; P2 f_gauge > f_livre nos 3 casos;
P3 f alinhado > desalinhado no gauge.  Resultado: P1 ✗ (2/3 casos), P2 ✗ (m=0), P3 ✗ (J=2).
Achados: (i) a IM de blocos do gauge também trava (~8 sítios; livre 14–18) e decai monotonicamente,
sem o pico balístico de pares de quasipartículas no cone de luz que o livre exibe;
(ii) a fração multipartida f do gauge cresce com o confinamento (0,57 → 0,75 → 0,86).
"""
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("tenpy")
from pesquisa_fisica.emergence.quench import _ptrace, _S_rho, block_pair_rdms  # noqa: E402

RES = Path(__file__).resolve().parents[2] / "results"
CASES = [(0.0, 1.0), (0.5, 1.0), (0.5, 2.0)]


def _R():
    return {(r["m"], r["J"]): r for r in json.load(open(RES / "multipartido_N24.json"))}


def _reach(rows, key, thr=1e-3):
    ds = [r["d"] for r in rows if r[key] > thr]
    return max(ds) if ds else 0


def _f(rows):
    v = [(q["I_block"] - q["I_pairs"]) / q["I_block"] for q in rows if q["d"] >= 4 and q["I_block"] > 1e-4]
    return float(np.mean(v))


def test_mdr_de_blocos_disjuntos_reproduz_im_de_correlatores():
    from pesquisa_fisica.emergence.schwinger_dmrg import gauge_profiles, ground_and_first_excited
    psi = ground_and_first_excited(14, 1.0, 0.3, 0.8, chi_max=40, excited=False)["psi0"]
    rh = block_pair_rdms(psi, 3, [6, 9], size=2)
    g = gauge_profiles(14, 0.8, 3, [6, 9], m=0.3, chi_max=40)
    for k, b in enumerate([6, 9]):
        R = rh[b]
        assert np.trace(R).real == pytest.approx(1.0, abs=1e-10)
        mi = _S_rho(_ptrace(R, [0], 4)) + _S_rho(_ptrace(R, [2], 4)) - _S_rho(_ptrace(R, [0, 2], 4))
        assert mi == pytest.approx(g["MI_s"][k], abs=1e-8)


def test_P1_falhou_ganho_de_alcance_nao_e_maior_no_gauge():
    R = _R()
    gain = lambda x: _reach(x["alinhado"], "I_block") - _reach(x["site"], "I")
    wins = [gain(R[c]["gauge"]) > gain(R[c]["free"]) for c in CASES]
    assert sum(wins) < len(CASES)                                 # falhou em pelo menos um caso (2 de 3)


def test_blocos_tambem_travam_no_gauge():
    R = _R()
    for c in CASES:
        assert _reach(R[c]["gauge"]["alinhado"], "I_block") <= 8
        assert _reach(R[c]["free"]["alinhado"], "I_block") >= 14


def test_gauge_sem_pico_balistico_livre_com_pico():
    """Livre: IM de blocos NÃO monótona (pares de quasipartículas no cone de luz). Gauge: decai monotonicamente."""
    R = _R()
    for c in CASES:
        g = np.array([q["I_block"] for q in R[c]["gauge"]["alinhado"]])
        f = np.array([q["I_block"] for q in R[c]["free"]["alinhado"]])
        assert np.all(np.diff(g[g > 1e-4]) < 0)
        assert np.any(np.diff(f[f > 1e-4]) > 0)


def test_P2_falhou_em_m0_mas_f_cresce_com_confinamento():
    R = _R()
    fg = [_f(R[c]["gauge"]["alinhado"]) for c in CASES]
    ff = [_f(R[c]["free"]["alinhado"]) for c in CASES]
    assert fg[0] < ff[0]                                          # m=0: P2 falha
    assert fg[1] > ff[1] and fg[2] > ff[2]
    assert fg[0] < fg[1] < fg[2]                                  # f cresce com o confinamento


def test_P3_falhou_em_J2():
    R = _R()
    al = [_f(R[c]["gauge"]["alinhado"]) for c in CASES]
    de = [_f(R[c]["gauge"]["desalinhado"]) for c in CASES]
    assert al[0] > de[0] and al[1] > de[1]
    assert al[2] <= de[2]
