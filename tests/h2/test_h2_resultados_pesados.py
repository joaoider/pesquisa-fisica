"""Verifica os resultados gravados pelos cálculos pesados (workflow "Cálculo pesado" no GitHub Actions).

Não recalcula nada: lê results/*.json e confere as predições pré-registradas. Se o arquivo não existir
(cálculo ainda não rodado), o teste é pulado.
"""
import json
from pathlib import Path

import numpy as np
import pytest

RESULTS = Path(__file__).resolve().parents[2] / "results"


def carregar(nome):
    p = RESULTS / nome
    if not p.exists():
        pytest.skip(f"{nome} ainda não foi gerado pelo workflow 'Cálculo pesado'")
    return json.loads(p.read_text())


def I_media_por_distancia(I):
    I = np.asarray(I)
    n = len(I)
    return np.array([np.mean([I[i, i + d] for i in range(n - d)]) for d in range(1, n)])


def test_su3_n8_predicoes_confirmadas():
    d = carregar("su3_n8.json")
    assert d["N"] == 8 and d["Nc"] == 3
    v = d["veredito"]
    assert v["P1_fermions_livres"] and v["P1_singleto"]
    assert v["P2_cadeia_para_todo_J"]
    assert v["P3_localidade_cresce_com_J"]


def test_su3_n8_sem_forca_zigue_zague_com_forca_monotonico():
    """Achado da v0.8 confirmado em N = 8: em J = 0 a informação mútua média oscila com a paridade da
    distância; com J > 0 ela decai monotonicamente."""
    d = carregar("su3_n8.json")
    livre = I_media_por_distancia(d["por_J"]["0.0"]["informacao_mutua"])
    assert livre[2] > livre[1]                      # I(d=3) > I(d=2)
    for J in ("0.5", "2.0", "5.0"):
        perfil = I_media_por_distancia(d["por_J"][J]["informacao_mutua"])
        assert np.all(np.diff(perfil[:4]) < 0)      # decai nas 4 primeiras distâncias
