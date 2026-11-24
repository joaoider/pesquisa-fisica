"""Ponte H1↔H2 (v0.29): Z_N (N=2..5) em 2+1D, rede 3×3 aberta, h=0,5, K=1 — fração destilável × N.
Pré-registro: P1 fração média monótona em N; P2 N=2 ≤ 0,05 e N=5 ≥ 0,15; P3 lei de perímetro (R² ≥ 0,9).
Resultado: P1 ✓ no conjunto pré-registrado (todas as regiões tratáveis), mas inversão N=4→5 no conjunto comum;
P2 ✗ (N=5: 0,001 — duas ordens abaixo); P3 ✓. N=4≈5: quase convergido ao U(1) em acoplamento intermediário
(⟨E²⟩≈0,37, P(|E|=N/2)≈3%); regime de fótons (acoplamento fraco) fora do alcance de diagonalização exata.
"""
import json
from pathlib import Path

import numpy as np
import pytest

from pesquisa_fisica.emergence.z2_links import Z2Lattice, region_entropies
from pesquisa_fisica.emergence.zn_links import ZNLattice, region_entropies_zn

RES = Path(__file__).resolve().parents[1] / "results"


def _D():
    return {d["summary"]["N"]: d for d in json.load(open(RES / "zn_fracao_bulk_h0.5.json"))}


def test_N2_reproduz_Z2():
    h2 = 0.3
    z = Z2Lattice(3, 3); Ez, _, pz = z.ground_state(h2)
    zn = ZNLattice(3, 3, 2); En, _, pn = zn.ground_state(4 * h2)
    assert En - h2 * z.L == pytest.approx(Ez, abs=1e-8)
    E = zn.link_values(zn.digits())
    for r in [(0, 1, 0, 1), (0, 2, 0, 2)]:
        a = region_entropies(z, pz, z.rect_links(*r)); b = region_entropies_zn(zn, pn, zn.rect_links(*r), E)
        assert a["S_ext"] == pytest.approx(b["S_ext"], abs=1e-8) and a["H_edge"] == pytest.approx(b["H_edge"], abs=1e-8)


def test_gauss_automatica():
    lat = ZNLattice(2, 3, 3)
    E = lat.link_values(lat.digits())
    _, inc = lat.cut_vertices([0])
    for v, ks in inc.items():
        d = sum((E[:, k] if lat.links[k][0] == v else -E[:, k]) for k in ks)
        assert np.all(d % 3 == 0)


def test_P1_monotona_no_conjunto_pre_registrado():
    D = _D()
    f = [D[N]["summary"]["frac_dist_mean"] for N in (2, 3, 4, 5)]
    assert np.all(np.diff(f) > 0)


def test_P2_falhou_borda_continua_dominando():
    D = _D()
    assert D[2]["summary"]["frac_dist_mean"] <= 0.05
    assert D[5]["summary"]["frac_dist_mean"] < 0.15


def test_P3_lei_de_perimetro():
    for d in _D().values():
        assert d["summary"]["R2_edge"] >= 0.9
