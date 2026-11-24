"""H2-T1 controle com interação local (v0.23): I₃_inv após o quench — gauge × H_V = hopping + massa + V Σ n n
(mesmo gap, exato N=20). Pré-registro: P1 controles com I₃_inv > −0,01 (≈0); P2 |I₃_gauge| ≥ 3·|I₃_controle|.
Resultado: P1 ✗ (V=1,2 no caso m=0: −0,015/−0,016); P2 ✓ (3,8–16× em t=6; 3,7–220× na média t∈[3,6]).
Exploratório: com confinamento, pares longos I_inv(L:R) do gauge 13–53× menores que os dos controles locais.
"""
import json
from pathlib import Path

import numpy as np

from pesquisa_fisica.emergence.ed_quench import evolve, evolve_H, local_interaction_hamiltonian

RES = Path(__file__).resolve().parents[2] / "results"


def _data():
    ref = {(r["m"], r["J"]): r for r in json.load(open(RES / "quench_invariante_N20.json"))}
    return ref, json.load(open(RES / "controle_interacao_N20.json"))


def test_controle_V0_igual_ao_livre():
    N = 12
    H, b = local_interaction_hamiltonian(N, 1.0, 0.5, 0.0)
    assert np.abs(evolve_H(H, b, N, [1.0])[0] - evolve(N, 0.5, 0.0, [1.0])[0]).max() < 1e-12


def test_gaps_casados():
    _, C = _data()
    for c in C:
        assert abs(c["gap"] - c["gap_target"]) < 1e-3


def test_P1_falhou_interacao_local_tambem_da_I3_negativa():
    _, C = _data()
    assert any(c["I3_i"][-1] <= -0.01 for c in C)


def test_P2_passou_gauge_muito_maior():
    ref, C = _data()
    for c in C:
        g = ref[tuple(c["gauge_case"])]["gauge"]["I3_i"][-1]
        assert abs(g) >= 3 * abs(c["I3_i"][-1])


def test_media_temporal_robusta():
    ref, C = _data()
    t = np.array(C[0]["times"]); w = t >= 3
    for c in C:
        g = np.array(ref[tuple(c["gauge_case"])]["gauge"]["I3_i"])[w].mean()
        assert abs(g) >= 3 * abs(np.array(c["I3_i"])[w].mean())


def test_exploratorio_confinamento_mata_pares_longos_alem_da_interacao_local():
    ref, C = _data()
    g = ref[(0.5, 2.0)]["gauge"]["ILR_i"][-1]
    for c in C:
        if tuple(c["gauge_case"]) == (0.5, 2.0):
            assert c["ILR_i"][-1] > 10 * g
