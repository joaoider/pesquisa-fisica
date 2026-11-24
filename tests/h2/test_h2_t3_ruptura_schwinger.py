"""H2-T3 ressonância de ruptura no Schwinger real (v0.26): Ī_inv(D) entre células, quench exato N=20.
Pré-registro: P1 pico interior em D ∈ [4,12]; P2 D_pico acompanha d* = 2·gap/J. Resultado: P1 ✗, P2 ✗ —
Ī(D) monótona decrescente em todos os casos e em TODOS os instantes. A ressonância da v0.25 era artefato do
corte (força nula além de R ⇒ partículas livres); a ruptura física gera mésons ligados, sem transporte longo.
"""
import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parents[2] / "results"
DS = [2, 4, 6, 8, 10, 12, 14]


def _rows():
    return json.load(open(RES / "ruptura_schwinger_N20.json"))


def test_casos_cobrem_d_star_de_3_a_9():
    d = sorted(r["d_star"] for r in _rows())
    assert d[0] < 4 and d[-1] > 8


def test_P1_falhou_sem_pico_interior():
    for r in _rows():
        I = np.array([r["Ibar"][str(D)] for D in DS])
        sig = I > 1e-5
        assert np.all(np.diff(I[sig]) < 0)


def test_monotona_em_todos_os_instantes():
    for r in _rows():
        M = np.array([r["I"][str(D)] for D in DS])
        for k in range(M.shape[1]):
            c = M[:, k][M[:, k] > 1e-4]
            assert np.all(np.diff(c) <= 0.1 * c[:-1] + 1e-12)
