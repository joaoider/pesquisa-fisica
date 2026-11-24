"""H2-T1 potencial linear truncado em R, sem gauge (v0.24). Identidade: J Σ L_n² = −(J/2) Σ Q_k Q_l |k−l| (Q_tot=0).
Pré-registro: P0 H_{N−1} = Schwinger; P1 |I₃| monótono em R; P2 I(L:R) confinante decrescente em R;
P3 R=8 ≥ 80% do |I₃| do gauge. Resultado: P0 ✓ (exato), P1 ✗ (caso confinante: R=5,8 passam do gauge),
P2 ✗ (só R = alcance total suprime os pares), P3 ✓.
Leitura: cada efeito é controlado pelo potencial até a ESCALA do observável — I₃ (6 sítios) satura com R ≳ 5–8;
a supressão de pares separados por 8 sítios só aparece com R > 8 (alcance completo, como a lei de Gauss impõe).
"""
import json
from pathlib import Path

import numpy as np

from pesquisa_fisica.emergence.ed_quench import truncated_linear_hamiltonian
from pesquisa_fisica.emergence.schwinger import hamiltonian

RES = Path(__file__).resolve().parents[2] / "results"


def _rows(case):
    return {r["R"]: r for r in json.load(open(RES / "alcance_linear_N20.json")) if tuple(r["gauge_case"]) == case}


def test_P0_identidade_gauss_potencial_linear():
    for N in (8, 12):
        for m, J in [(0.0, 1.0), (0.5, 2.0), (0.3, 0.7)]:
            assert abs(truncated_linear_hamiltonian(N, 1.0, m, J, N - 1)[0] - hamiltonian(N, w=1.0, m=m, J=J)[0]).max() < 1e-10


def test_gaps_casados():
    for case in [(0.0, 1.0), (0.5, 2.0)]:
        for r in _rows(case).values():
            assert abs(r["gap"] - r["gap_target"]) < 1e-3


def test_P1_monotono_sem_massa_mas_falhou_no_confinante():
    a = [abs(_rows((0.0, 1.0))[R]["I3_avg"]) for R in (2, 3, 5, 8, 19)]
    assert np.all(np.diff(a) > 0)
    b = [abs(_rows((0.5, 2.0))[R]["I3_avg"]) for R in (2, 3, 5, 8, 19)]
    assert not np.all(np.diff(b) > 0)                          # falha pré-registrada
    assert b[3] > b[4]                                         # R=8 passa do gauge


def test_P2_falhou_supressao_de_pares_so_com_alcance_total():
    c = [_rows((0.5, 2.0))[R]["ILR_i"][-1] for R in (2, 3, 5, 8, 19)]
    assert not np.all(np.diff(c) < 0)
    assert c[-1] < 0.05 * min(c[:-1])                          # só R completo suprime (≥ 20×)
    d = [_rows((0.0, 1.0))[R]["ILR_i"][-1] for R in (2, 3, 5, 8, 19)]
    assert d[-1] < 0.5 * min(d[:-1])


def test_P3_R8_reproduz_I3_do_gauge():
    for case in [(0.0, 1.0), (0.5, 2.0)]:
        r = _rows(case)
        assert abs(r[8]["I3_avg"]) >= 0.8 * abs(r[19]["I3_avg"])
