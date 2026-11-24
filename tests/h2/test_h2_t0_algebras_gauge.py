"""H2-T0 localidade em álgebras invariantes de gauge (v0.20) — Casini-Huerta-Rosabal em Schwinger 1+1D.

Álgebras: qubit (antiga) | inv (comuta com Q de cada região; sem Wilson entre regiões) | ce (inv + fluxos de
fronteira) | carga (diagonal).  Pré-registro: P1 inv causal (t_pico crescente, v_b ∈ [1,5;2,3]);
P2 v_ce > v_inv + 0,5; P3 razão I(d=1)/I(d=3) do vácuo cresce com J em inv.
Resultado: P1 ✓ (v=1,89), P2 ✗ (v_ce=1,76: o centro elétrico também é causal), P3 ✓.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from pesquisa_fisica.emergence.ed_quench import block_mi, classical_block_mi, evolve, symmetric_blocks
from pesquisa_fisica.emergence.gauge_algebra import algebra_entropy, algebra_mi

RES = Path(__file__).resolve().parents[2] / "results"
RS = (2, 4, 6, 8)


def _d():
    return json.load(open(RES / "algebras_gauge.json"))


def test_modos_consistentes_com_funcoes_anteriores():
    N = 12
    s = evolve(N, 0.0, 1.0, [1.5])[0]
    A, B = symmetric_blocks(N, 2, 2)
    assert algebra_mi(s, N, A, B, "qubit") == pytest.approx(block_mi(s, N, A, B), abs=1e-12)
    assert algebra_mi(s, N, A, B, "carga") == pytest.approx(classical_block_mi(s, N, A, B), abs=1e-12)
    assert algebra_entropy(s, N, [A], [], "ce") == pytest.approx(algebra_entropy(s, N, [A], [], "inv"), abs=1e-12)
    s0 = evolve(N, 0.0, 1.0, [0.0])[0]
    for md in ("qubit", "inv", "ce", "carga"):
        assert abs(algebra_mi(s0, N, A, B, md)) < 1e-12


def test_hierarquia_carga_inv():
    """inv ⊇ carga ⇒ I_inv ≥ I_carga (monotonicidade da IM sob restrição de álgebra)."""
    q = _d()["quench"]
    for r in RS:
        assert np.all(np.array(q[f"inv_{r}"]["I"]) >= np.array(q[f"carga_{r}"]["I"]) - 1e-10)


def test_P1_passou_inv_e_causal():
    q = _d()["quench"]
    tp = [q[f"inv_{r}"]["t_peak"] for r in RS]
    assert np.all(np.diff(tp) > 0)
    assert 1.5 <= q["inv_front"]["v"] <= 2.3
    assert q["qubit_front"]["v"] > 3.0                       # contraste: a álgebra de qubits não é causal


def test_P2_falhou_centro_eletrico_tambem_causal():
    q = _d()["quench"]
    assert q["ce_front"]["v"] <= q["inv_front"]["v"] + 0.5
    assert 1.5 <= q["ce_front"]["v"] <= 2.3


def test_P3_passou_localidade_cresce_com_J_em_inv():
    v = _d()["vacuo"]
    r = [v[J]["inv"]["ratio"] for J in ("0.0", "0.5", "2.0")]   # J=5: I(d=3) ~ 1e-11, no piso numérico
    assert r[0] < r[1] < r[2]


def test_baseline_J0_muito_mais_local_na_algebra_invariante():
    v = _d()["vacuo"]["0.0"]
    assert v["qubit"]["ratio"] < 100 and v["inv"]["ratio"] > 1000
