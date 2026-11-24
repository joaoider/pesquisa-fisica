"""H2-T1 razão gauge/livre (mesmo gap) na álgebra invariante (v0.21): vácuo Schwinger m=0 vs livre, DMRG N=64,
IM entre células de 2 sítios. Pré-registro: P1 razão inv ∈ [0,85;1,15]; P2 κ/gap ∈ [0,8;1,2]; P3 controle de
qubits ∈ [1,3;1,7]. Resultado: P1 ✓ (1,03–1,07), P2 ✓, P3 ✓ (1,46–1,48).
Conclusão: o fator ≈1,5 da v0.14–v0.15 era artefato de correlações carregadas; a mesmo gap, gauge = livre.
"""
import json
from pathlib import Path

import numpy as np
import pytest

pytest.importorskip("tenpy")
from pesquisa_fisica.emergence.gauge_algebra import algebra_mi, two_cell_mi  # noqa: E402
from pesquisa_fisica.emergence.quench import block_pair_rdms  # noqa: E402
from pesquisa_fisica.emergence.schwinger import ground_state  # noqa: E402
from pesquisa_fisica.emergence.schwinger_dmrg import ground_and_first_excited  # noqa: E402

RES = Path(__file__).resolve().parents[2] / "results"


def _rows():
    return json.load(open(RES / "razao_invariante_N64.json"))


def test_two_cell_mi_mps_igual_a_diagonalizacao_exata():
    N, m, J = 12, 0.5, 1.0
    psi = ground_and_first_excited(N, 1.0, m, J, chi_max=60, excited=False)["psi0"]
    _, full, _ = ground_state(N, w=1.0, m=m, J=J)
    for b in (6, 8):
        d = two_cell_mi(block_pair_rdms(psi, 2, [b], size=2)[b], (2, b))
        for md in ("qubit", "inv", "carga"):
            assert d[md] == pytest.approx(algebra_mi(full, N, [2, 3], [b, b + 1], md), abs=1e-8)


def test_P1_razao_invariante_igual_a_1():
    for r in _rows():
        assert 0.85 <= r["k_inv"]["ratio"] <= 1.15, r["J"]


def test_P2_taxa_inv_igual_ao_gap():
    for r in _rows():
        for side in ("g_over_gap", "f_over_gap"):
            assert 0.8 <= r["k_inv"][side] <= 1.2


def test_P3_controle_qubits_reproduz_v014():
    for r in _rows():
        assert 1.3 <= r["k_qubit"]["ratio"] <= 1.7


def test_inv_dominado_por_cargas():
    for r in _rows():
        assert abs(r["k_inv"]["ratio"] - r["k_carga"]["ratio"]) < 0.02
