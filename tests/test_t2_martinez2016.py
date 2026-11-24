"""T2 (validação da ferramenta): nosso Schwinger contra o experimento de Martinez et al., Nature 534, 516 (2016).
N=4 íons, J=w, m=0,5w, ε₀=0, quench do vácuo nu; dados digitalizados das figuras vetoriais (SVG) do ar5iv.
Pré-registro (previsão commitada antes de ler as figuras, 3339280): P1 curva ideal = a do artigo (|Δν| ≤ 0,03);
P2 dados a |Δν| ≤ 0,07 do IDEAL; P3 E_N medido > 0 e ≤ ideal + 0,1.
Resultado: P1 ✓ (0,0003), P2 ✗ (+0,077 e +0,173 nos dois últimos passos), P3 ✓.
Exploratório: Trotter Z·pares·ímpares (4 passos de π/8) + desfase p=0,038 reproduz EXATAMENTE os pontos de
discretização e o modelo de erro do artigo; rms contra os dados cai de 0,105 (ideal) para 0,015.
"""
import json
from pathlib import Path

import numpy as np

RES = Path(__file__).resolve().parents[1] / "results"


def _j(n):
    return json.load(open(RES / n))


def test_P1_curva_ideal_identica_a_do_artigo():
    a = _j("t2_martinez2016_avaliacao.json")
    assert a["P1"] and a["P1_max_dnu"] < 0.005


def test_P2_falhou_contra_o_ideal():
    a = _j("t2_martinez2016_avaliacao.json")
    d = np.abs(a["dnu_exp_menos_ideal"])
    assert not a["P2"] and d[3] > 0.07 and d[4] > 0.07
    assert d[1] < 0.02 and d[2] < 0.02                        # nos primeiros passos o ideal basta


def test_P3_emaranhamento_medido_abaixo_do_ideal():
    assert _j("t2_martinez2016_avaliacao.json")["P3"]


def test_exploratorio_reproduz_modelo_de_erro_do_artigo():
    x = _j("t2_martinez2016_exploratorio.json")
    D = _j("t2_martinez2016_dados_digitalizados.json")
    tz = x["esquemas"]["Trotter Z·par·ímpar"]
    assert np.max(np.abs(np.array(tz["sem_ruido"]) - x["red_artigo"])) < 1e-3
    assert np.max(np.abs(np.array(tz["com_desfase"][1:]) - np.array(D["nu_modelo_erro_artigo"][1:]))) < 1e-3
    rms_ideal = np.sqrt(np.mean((np.array(x["esquemas"]["ideal (exato)"]["sem_ruido"]) - D["nu_exp"]) ** 2))
    rms_mod = np.sqrt(np.mean((np.array(tz["com_desfase"]) - D["nu_exp"]) ** 2))
    assert rms_mod < 0.02 and rms_ideal > 0.09
