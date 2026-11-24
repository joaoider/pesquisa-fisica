"""Ponte H1↔H2 (v0.28): centro elétrico em 2+1D — Z₂ pura, rede aberta 4×4 (40 elos), Hilbert estendido.
S_ext = H_borda (Shannon dos fluxos nos vértices de corte) + S_dest. Pré-registro: P1 desconfinado H_borda/S_ext ≥ 0,9;
P2 H_borda linear em n_corte (R² ≥ 0,95), inclinação decrescente em h; P3 intercepto ≈ −ln2 (desconf.) e ≈ 0 (conf.).
Resultado: P1 ✓ (0,998), P2 ✓ (0,983/0,953/0,952; 0,60 > 0,39 > 0,11), P3 ✓ (−0,48 ln2; −0,06 ln2).
"""
import json
from pathlib import Path

import numpy as np
import pytest

from pesquisa_fisica.emergence.z2_links import Z2Lattice, region_entropies

RES = Path(__file__).resolve().parents[1] / "results"


def _S():
    return {d["summary"]["h"]: d for d in json.load(open(RES / "centro_eletrico_2d.json"))}


def test_limite_toric_code_exato():
    """h → 0: S_ext = (n_corte − 1) ln 2, toda no centro elétrico (Hamma-Ionicioiu-Zanardi; Kitaev-Preskill)."""
    lat = Z2Lattice(3, 3)
    _, _, psi = lat.ground_state(1e-9)
    for r in [(0, 1, 0, 1), (1, 2, 1, 2), (0, 2, 0, 2), (1, 3, 0, 3)]:
        d = region_entropies(lat, psi, lat.rect_links(*r))
        assert d["S_ext"] == pytest.approx((d["n_cut"] - 1) * np.log(2), abs=1e-8)
        assert abs(d["S_dist"]) < 1e-8


def test_S_dest_nao_negativa():
    for d in _S().values():
        assert d["summary"]["S_dist_min"] > -1e-10


def test_P1_desconfinado_centro_carrega_quase_tudo():
    assert _S()[0.15]["summary"]["edge_frac_min"] >= 0.9


def test_P2_lei_de_perimetro_do_centro():
    s = [_S()[h]["summary"] for h in (0.15, 0.33, 0.8)]
    assert all(x["R2_H"] >= 0.95 for x in s)
    assert s[0]["slope_H"] > s[1]["slope_H"] > s[2]["slope_H"]


def test_P3_entropia_topologica():
    S = _S()
    assert -1.3 <= S[0.15]["summary"]["b_S_over_ln2"] <= -0.4
    assert abs(S[0.8]["summary"]["b_S_over_ln2"]) <= 0.3


def test_achado_centro_domina_em_todas_as_fases():
    for d in _S().values():
        assert d["summary"]["edge_frac_min"] > 0.97
