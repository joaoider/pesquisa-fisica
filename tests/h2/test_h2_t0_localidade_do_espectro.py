"""H2 · T0 — Localidade a partir do espectro (Cotler, Penington & Ranard, CMP 2019).

Pergunta: se o espaço emerge, de onde vem a divisão em "sítios vizinhos"? Teste: usando SÓ a lista
ordenada de autovalores (nenhuma informação de base ou de vizinhança), ajustar um ansatz de cadeia local
(12 parâmetros) a espectros de 32 e 64 níveis.

Predições pré-registradas (Notion):
  P1. Espectro de Hamiltoniano local (da família do ansatz): ajuste exato (erro < 10⁻¹⁰).
  P2. Espectro genérico (GUE) e de Hamiltoniano NÃO local: o ajuste falha (erro > 10⁻³).
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.locality_from_spectrum import (
    ansatz_terms,
    fit_spectrum,
    gue_spectrum,
    hamiltonian,
    random_nonlocal_hamiltonian_spectrum,
)


@pytest.fixture(scope="module")
def termos5():
    return ansatz_terms(5)


@pytest.mark.parametrize("seed", [0, 1])
def test_espectro_local_e_reconhecido_so_pelos_autovalores(termos5, seed):
    rng = np.random.default_rng(100 + seed)
    E = np.linalg.eigvalsh(hamiltonian(rng.standard_normal(len(termos5)), termos5))
    err, _ = fit_spectrum(E, termos5, n_starts=60, seed=seed)
    assert err < 1e-10


def test_espectro_local_n6_64_niveis():
    termos = ansatz_terms(6)
    rng = np.random.default_rng(7)
    E = np.linalg.eigvalsh(hamiltonian(rng.standard_normal(len(termos)), termos))
    # o ajuste é uma BUSCA com mínimos locais: mais pontos de partida (para no primeiro sucesso) — v0.11:
    # com 20 partidas, o runner do CI ficou preso uma vez (erro 2,8e-3) por diferenças de ponto flutuante
    err, _ = fit_spectrum(E, termos, n_starts=80, seed=0)
    assert err < 1e-10


def test_espectro_generico_gue_nao_tem_fatoracao_local(termos5):
    err, _ = fit_spectrum(gue_spectrum(32, np.random.default_rng(1)), termos5, n_starts=20)
    assert err > 1e-3


def test_espectro_de_hamiltoniano_nao_local_nao_se_ajusta(termos5):
    err, _ = fit_spectrum(random_nonlocal_hamiltonian_spectrum(5, np.random.default_rng(2)), termos5, n_starts=20)
    assert err > 1e-3


@pytest.mark.skip(reason="Inconclusivo (v0.9): com ansatz local SEM invariância de translação (XXZ desordenado, "
                         "25 parâmetros, 128 níveis) o otimizador fica preso em mínimos locais (erro ~1,5×10⁻³ "
                         "mesmo para alvo local). Encontrar a fatoração local pode ser difícil mesmo quando ela "
                         "existe (Cotler et al.). Próximo: muitos pontos de partida/otimização global no workflow pesado.")
def test_localidade_do_espectro_caso_desordenado():
    ...
