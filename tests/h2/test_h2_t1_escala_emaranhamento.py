"""H2-T1 escala do emaranhamento após o quench (v0.18): N=24, t=6, gauge vs livre de mesmo gap (TDVP, spin).

Pré-registro: P1 ℓ*_gauge < ℓ*_livre (saturação de S(ℓ)); P2 I(L:R)_gauge/I_livre < 0,3 em ℓ=8;
P3 I₃ (blocos de 2) gauge < 0 e < livre — nos 3 casos (m,J) = (0,1), (0,5;1), (0,5;2).
Resultado: P1 ✗ (m=0: empate ℓ*=9); P2 ✗ (m=0: 0,40); P3 ✗ (m=0: −0,107 > −0,161). Com m>0 as três passam.
Leitura: confinamento ⇒ emaranhamento de curto alcance (lei de área); blindagem (m=0) mantém pares longos
(I(L:R) em ℓ=8 é 40% do livre), invisíveis às regiões pequenas da v0.17.
"""
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("tenpy")
from pesquisa_fisica.emergence.quench import block_entropy, cut_entropy, tripartite_info  # noqa: E402

RES = Path(__file__).resolve().parents[2] / "results"
CASES = [(0.0, 1.0), (0.5, 1.0), (0.5, 2.0)]


def _R():
    return {(r["m"], r["J"]): r for r in json.load(open(RES / "escala_N24.json"))}


def _ell_star(x, frac=0.9):
    S = np.array([r["S_block"] for r in x["S_ell"]])
    return int(np.argmax(S >= frac * S[-1])) + 1


def _ILR(x, ell=8):
    return [r["I_LR"] for r in x["S_ell"] if r["ell"] == ell][0]


def test_entropia_de_bloco_confere_com_tenpy():
    from pesquisa_fisica.emergence.schwinger_dmrg import ground_and_first_excited
    psi = ground_and_first_excited(14, 1.0, 0.3, 0.8, chi_max=40, excited=False)["psi0"]
    assert block_entropy(psi, 0, 5) == pytest.approx(cut_entropy(psi, 5), abs=1e-9)
    assert block_entropy(psi, 3, 6) == pytest.approx(psi.entanglement_entropy_segment(list(range(3, 9)))[0], abs=1e-9)
    assert np.isfinite(tripartite_info(psi, 4, 2))


def test_P1_falhou_em_m0_empate():
    R = _R()
    assert _ell_star(R[(0.0, 1.0)]["gauge"]) == _ell_star(R[(0.0, 1.0)]["free"]) == 9   # falha pré-registrada
    for c in [(0.5, 1.0), (0.5, 2.0)]:
        assert _ell_star(R[c]["gauge"]) < _ell_star(R[c]["free"]), c


def test_P2_falhou_em_m0_pares_longos_sobrevivem_a_blindagem():
    R = _R()
    ratio = {c: _ILR(R[c]["gauge"]) / _ILR(R[c]["free"]) for c in CASES}
    assert ratio[(0.0, 1.0)] >= 0.3                          # falha pré-registrada
    assert ratio[(0.5, 1.0)] < 0.3 and ratio[(0.5, 2.0)] < 0.3
    assert ratio[(0.5, 2.0)] < ratio[(0.5, 1.0)] < ratio[(0.0, 1.0)]   # confinamento mata os pares longos


def test_P3_falhou_em_m0():
    R = _R()
    g = {c: R[c]["gauge"]["I3_s2"] for c in CASES}
    f = {c: R[c]["free"]["I3_s2"] for c in CASES}
    assert g[(0.0, 1.0)] > f[(0.0, 1.0)]                      # falha pré-registrada
    for c in [(0.5, 1.0), (0.5, 2.0)]:
        assert g[c] < 0 and g[c] < f[c]


def test_confinamento_produz_lei_de_area():
    """Com m>0 e J=2, S(ℓ) do gauge para de crescer (lei de área); o livre segue crescendo."""
    x = _R()[(0.5, 2.0)]
    Sg = np.array([r["S_block"] for r in x["gauge"]["S_ell"]])
    Sf = np.array([r["S_block"] for r in x["free"]["S_ell"]])
    assert Sg[-1] - Sg[5] < 0.05 and Sf[-1] - Sf[5] > 0.5
