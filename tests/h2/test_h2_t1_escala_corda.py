"""H2-T1 lei de escala da supressão de pares (v0.25): potencial linear truncado em R, gap casado (J=2, gap 3,45),
quench exato N=20, blocos de 4 sítios a distância centro-a-centro D ∈ {8,10,12,14}, Ī = média t∈[3,6] de I_inv.
Pré-registro: P1 colapso de f em R/D (degrau); P2 R_meio ∝ D. Resultado: P1 ✗, P2 ✗ (f não monótona).
Exploratório: RESSONÂNCIA — com R truncado a IM entre blocos é AMPLIFICADA (até 7× o caso local) e, para R ≥ 8,
o pico em D acompanha R (8→8, 10→10, 12→12); o gauge (R sem corte) suprime tudo. Efeito par/ímpar em R.
"""
import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parents[2] / "results"
RS = [4, 6, 8, 10]
DS = [r + 4 for r in RS]


def _by():
    return {r["R"]: r for r in json.load(open(RES / "escala_corda_N20.json"))}


def _f(by, R):
    g = {r: by[19]["Ibar"][str(r)] for r in RS}
    l = {r: by[2]["Ibar"][str(r)] for r in RS}
    return [(by[R]["Ibar"][str(r)] - g[r]) / (l[r] - g[r]) for r in RS]


def test_gaps_casados():
    for r in _by().values():
        assert abs(r["gap"] - r["gap_target"]) < 1e-3


def test_P1_falhou_sem_colapso_em_degrau():
    by = _by()
    viol = 0
    for R in [3, 4, 5, 6, 8, 10, 12]:
        for f, D in zip(_f(by, R), DS):
            x = R / D
            viol += (x <= 0.6 and f < 0.5) or (x >= 1.5 and f > 0.2)
    assert viol > 0


def test_P2_falhou_f_nao_monotona_em_R():
    by = _by()
    for k in range(len(RS)):
        col = [_f(by, R)[k] for R in [2, 3, 4, 5, 6, 8, 10, 12]]
        assert not np.all(np.diff(col) <= 0)


def test_ressonancia_pico_acompanha_R():
    by = _by()
    for R in (8, 10, 12):
        I = [by[R]["Ibar"][str(r)] for r in RS]
        assert DS[int(np.argmax(I))] == R
    assert max(_f(by, 12)) > 5                                  # amplificação forte


def test_gauge_suprime_em_todas_as_distancias():
    by = _by()
    for r in RS:
        assert by[19]["Ibar"][str(r)] < min(by[R]["Ibar"][str(r)] for R in [2, 3, 4, 5, 6, 8, 10, 12])
