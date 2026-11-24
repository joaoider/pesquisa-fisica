"""H2-T3 (v0.31): I₃ de Rényi-2 invariante como predição mensurável — gauge (m=0,J=1) × controle local (V=1, mesmo gap),
N=16, Trotter Z·par·ímpar dt=0,25 (ordem validada na v0.30), desfase p, sombras clássicas em lotes.
Pré-registro: P1 razão exata ≥ 3; P2 retenção de Δ ≥ 50% com p=0,01; P3 |Δ| ≥ 3σ com M=10⁴.
Resultado: P1 ✗ (2,53), P2 ✗ (42%), P3 ✗ (|Δ|/σ = 0,06; seriam ~10⁷ disparos). Assinatura NÃO é testável na prática.
"""
import json
from pathlib import Path

import numpy as np

from pesquisa_fisica.emergence.ed_quench import evolve
from pesquisa_fisica.emergence.protocolo_t3 import (I3_renyi2_inv_from_rhoABC, batch_purity_fn, diag_energies, rho_regions,
                                                    shadow_batches, trotter_trajectory)

RES = Path(__file__).resolve().parents[2] / "results"


def _d():
    return json.load(open(RES / "t3_protocolo_i3.json"))


def test_trotter_converge_e_ordem_validada():
    rng = np.random.default_rng(0)
    N = 8
    psi = trotter_trajectory(N, diag_energies(N, 0.5, J=1.0), 0.01, 200, 0.0, rng)
    assert abs(np.vdot(psi, evolve(N, 0.5, 1.0, [2.0])[0])) ** 2 > 0.999
    st = trotter_trajectory(4, diag_energies(4, 0.5, J=1.0), np.pi / 8, 4, 0.0, rng, record={1, 2, 3, 4})
    idx = np.arange(16); b = (idx[:, None] >> np.arange(4)) & 1
    exc = np.where(np.arange(4) % 2 == 0, b, 1 - b).sum(1) / 4
    nu = [float((np.abs(st[k]) ** 2 * exc).sum()) for k in (1, 2, 3, 4)]
    assert np.allclose(nu, [0.1982, 0.4716, 0.4235, 0.2252], atol=1e-3)     # = Martinez 2016 (discretização)


def test_sombras_sem_vies():
    rng = np.random.default_rng(1)
    N = 8
    regs = {"A": [2, 3], "B": [4, 5], "C": [6, 7]}
    rho = rho_regions(evolve(N, 0.0, 1.0, [3.0])[0], N, regs)
    ex = I3_renyi2_inv_from_rhoABC(rho, regs)
    est = I3_renyi2_inv_from_rhoABC(rho, regs, batch_purity_fn(shadow_batches(rho, 6, 20000, 10, rng), 6))
    assert abs(est - ex) < 0.15


def test_gaps_casados():
    d = _d()
    assert abs(d["gap"] - d["gap_controle"]) < 1e-3


def test_P1_falhou_razao_abaixo_de_3():
    p = _d()["P1"]
    assert p["avg_gauge"] < 0 and abs(p["avg_gauge"]) < 3 * abs(p["avg_controle"])


def test_P2_falhou_retencao_abaixo_de_50():
    assert _d()["P2"]["retencao"] < 0.5


def test_P3_falhou_ruido_estatistico_domina():
    p = _d()["P3"]
    assert p["z"] < 3 and p["M_necessario_3sigma"] > 1e6
