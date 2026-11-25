"""H2-T3 triagem de observáveis baratos (v0.32): gauge × controles locais V∈{1,2,4} (mesmo gap), N=16.
Pré-registro: P1 F_gauge ≤ 0,5·min F_controle; P2 S₂(bloco 4) do gauge a ±30% de algum controle; P3 separação de F ≥ 3σ
com protocolo (Trotter dt=0,25, p=0,01, K=10⁴). Resultado: P1 ✗ (V=4 suprime mais), P2 ✓, P3 ✓ (11σ–40σ, mas sem
assinatura, pois P1 falhou). Exploratório: ν fica FORA da faixa de todos os controles nos dois casos.
"""
import json
from pathlib import Path

RES = Path(__file__).resolve().parents[2] / "results"


def _C():
    return {(c["m"], c["J"]): c["models"] for c in json.load(open(RES / "triagem_baratos_N16.json"))}


CTRL = ("V=1", "V=2", "V=4")


def test_P1_falhou_repulsao_local_forte_suprime_mais():
    for M in _C().values():
        assert M["gauge"]["exato"]["F"] > 0.5 * min(M[c]["exato"]["F"] for c in CTRL)


def test_P2_entropia_barata_nao_distingue():
    for M in _C().values():
        g = M["gauge"]["exato"]["S2_blk4"]
        assert any(abs(g - M[c]["exato"]["S2_blk4"]) <= 0.3 * g for c in CTRL)


def test_P3_F_mensuravel():
    for M in _C().values():
        g = M["gauge"]
        c = min(CTRL, key=lambda k: abs(M[k]["F_protocolo"] - g["F_protocolo"]))
        sig = (g["F_est_sigma"] ** 2 + M[c]["F_est_sigma"] ** 2) ** 0.5
        assert abs(g["F_protocolo"] - M[c]["F_protocolo"]) >= 3 * sig


def test_exploratorio_nu_fora_da_faixa_dos_controles():
    for M in _C().values():
        assert M["gauge"]["exato"]["nu"] > max(M[c]["exato"]["nu"] for c in CTRL)
