"""H2-T1 2+1D (v0.27): gauge Z₂ pura via dual TFIM (5×4 plaquetes), IM entre plaquetes nas álgebras full × inv.
Pré-registro: P1 confinado: full longe/perto ≥ 0,1 e inv ≤ 1e-3; P2 desconfinado: ambas ≤ 1e-2;
P3 Spearman(inv) ≥ 0,8 em todo h e Spearman(full, confinado) < 0,5.
Resultado: P1 ✓ (0,69 × 3e-7), P2 ✓, P3 ✗ na margem (inv ✓ 0,94–0,98; full confinado = 0,574 ≥ 0,5).
"""
import json
from functools import reduce
from pathlib import Path

import numpy as np

from pesquisa_fisica.emergence.z2_gauge_2d import bonds, ground_state

RES = Path(__file__).resolve().parents[2] / "results"


def _R():
    return {r["h"]: r for r in json.load(open(RES / "z2_2d_algebras.json"))}


def test_dual_esparso_igual_denso():
    X = np.array([[0, 1], [1, 0.]]); Z = np.diag([1., -1]); I = np.eye(2)
    op = lambda o, i, n: reduce(np.kron, [o if k == i else I for k in range(n)])
    Lx, Ly, h = 2, 3, 0.5
    n = Lx * Ly
    H = sum(-op(X, i, n) for i in range(n)) + sum(-h * op(Z, i, n) @ op(Z, j, n) for i, j in bonds(Lx, Ly))
    P = reduce(lambda a, b: a @ b, [op(X, i, n) for i in range(n)])
    w, V = np.linalg.eigh(H)
    ev = [w[k] for k in range(len(w)) if abs(V[:, k] @ P @ V[:, k] - 1) < 1e-6]
    assert abs(min(ev) - ground_state(Lx, Ly, h)[0]) < 1e-10


def test_h0_produto():
    E, gap, _ = ground_state(3, 3, 0.0)
    assert abs(E + 9) < 1e-10 and abs(gap - 4) < 1e-8


def test_P1_confinado_cordas_eletricas_de_longo_alcance():
    r = _R()[0.8]
    assert r["full"]["far_over_near"] >= 0.1
    assert r["inv"]["far_over_near"] <= 1e-3


def test_P2_desconfinado_ambas_decaem():
    r = _R()[0.15]
    assert r["full"]["far_over_near"] <= 1e-2 and abs(r["inv"]["far_over_near"]) <= 1e-2


def test_P3_inv_reconstroi_2D_em_todas_as_fases_full_degrada_mas_nao_abaixo_de_05():
    R = _R()
    for h in (0.15, 0.33, 0.8):
        assert R[h]["inv"]["spearman"] >= 0.8
    assert R[0.8]["full"]["spearman"] >= 0.5                       # falha pré-registrada (0,574)
    assert R[0.8]["full"]["spearman"] < R[0.8]["inv"]["spearman"] - 0.3
