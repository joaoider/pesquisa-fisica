"""H2-T1 confinamento na dinâmica na álgebra invariante (v0.22): quench exato N=20.

Pré-registro: P0 S de região única idêntica (qubit = inv); P1 m=0: razão gauge/livre de I(L:R) maior em inv
que em qubit; P2 m=0,5,J=2: razão inv < 0,3; P3 m=0,5,J=2: I₃_inv(gauge) < 0 e < livre.  Tudo ✓.
Exploratório: em inv, I₃ do livre ≈ 0 (quadro de pares) e a do gauge é negativa também em m=0.
"""
import json
from pathlib import Path

import numpy as np

from pesquisa_fisica.emergence.ed_quench import evolve
from pesquisa_fisica.emergence.gauge_algebra import algebra_entropy, algebra_I3

RES = Path(__file__).resolve().parents[2] / "results"


def _R():
    return {(r["m"], r["J"]): r for r in json.load(open(RES / "quench_invariante_N20.json"))}


def test_P0_teorema_regiao_unica():
    """|ψ⟩ com carga total definida ⇒ ρ_R comuta com Q_R ⇒ S_inv(R) = S_qubit(R)."""
    N = 12
    s = evolve(N, 0.3, 1.0, [2.0])[0]
    for R in ([0, 1, 2, 3, 4, 5], [3, 4, 5, 6, 7], [5]):
        assert abs(algebra_entropy(s, N, [R], (), "qubit") - algebra_entropy(s, N, [R], (), "inv")) < 1e-10
    for r in _R().values():
        for side in ("gauge", "free"):
            x = r[side]
            assert np.max(np.abs(np.array(x["S_half_q"]) - x["S_half_i"])) < 1e-10
            assert np.max(np.abs(np.array(x["S_M_q"]) - x["S_M_i"])) < 1e-10


def _ratio(r, k):
    return r["gauge"][k][-1] / r["free"][k][-1]


def test_P1_mesons_neutros_sobrevivem_a_algebra_invariante():
    r = _R()[(0.0, 1.0)]
    assert _ratio(r, "ILR_i") > _ratio(r, "ILR_q")


def test_P2_confinamento_mata_pares_longos_tambem_em_inv():
    assert _ratio(_R()[(0.5, 2.0)], "ILR_i") < 0.3


def test_P3_I3_invariante_negativa_e_abaixo_do_livre_com_confinamento():
    r = _R()[(0.5, 2.0)]
    assert r["gauge"]["I3_i"][-1] < 0 and r["gauge"]["I3_i"][-1] < r["free"]["I3_i"][-1]


def test_exploratorio_I3_livre_nula_e_gauge_negativa():
    for r in _R().values():
        assert abs(r["free"]["I3_i"][-1]) < 0.01
        assert r["gauge"]["I3_i"][-1] < -0.03


def test_I3_invariante_bem_definida():
    N = 10
    s = evolve(N, 0.0, 1.0, [0.0])[0]
    assert abs(algebra_I3(s, N, [3], [4], [5], "inv")) < 1e-12
