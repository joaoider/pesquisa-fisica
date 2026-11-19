"""H2 · T1 — Férmions em 2D num fundo de gauge U(1) (campo médio da QED₃), exato em rede 20×20 (400 sítios).

Predições pré-registradas (Notion):
  P1. Invariância de gauge U(1): uma transformação c_j → e^{iθ_j} c_j aleatória não muda o fluxo nem a
      informação mútua.
  P2. Em todos os fundos de fluxo (Φ = 0 metal, Φ = π Dirac, Hofstadter Φ = 2π/4 e 2π/5), o grafo dos
      4 parceiros mais emaranhados é a rede quadrada real.
  P3. O fluxo muda a INTENSIDADE do emaranhamento (a escala da métrica), não a vizinhança.
  P4. Fluxo aleatório (desordem de gauge) não destrói a geometria.
  P5. Controle: grafo aleatório não tem geometria.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.geometry_from_entanglement import graph_distances, knn_graph, spearman, upper
from pesquisa_fisica.emergence.u1_background import (
    gauge_transform_links,
    ground_state_correlations,
    hopping_matrix,
    landau_gauge_links,
    plaquette_flux,
    random_flux_links,
    site_mutual_information,
    torus_distances,
)

L = 20
T = torus_distances(L)
REDE = np.isclose(T, 1.0)
FLUXOS = {"metal_0": 0.0, "dirac_pi": np.pi, "hofstadter_1/4": 2 * np.pi / 4, "hofstadter_1/5": 2 * np.pi / 5}


def resolver(links):
    C, gap = ground_state_correlations(hopping_matrix(L, links))
    return C, site_mutual_information(C)


@pytest.fixture(scope="module")
def fundos():
    return {nome: resolver(landau_gauge_links(L, phi, twist=0.3)) for nome, phi in FLUXOS.items()}


# ---------------------------------------------------------------- P1

def test_invariancia_de_gauge_u1():
    links = landau_gauge_links(L, np.pi, twist=0.3)
    theta = np.random.default_rng(1).uniform(0, 2 * np.pi, L * L)
    links2 = gauge_transform_links(links, theta)
    assert np.allclose(np.exp(1j * plaquette_flux(L, links)), np.exp(1j * plaquette_flux(L, links2)))
    _, I1 = resolver(links)
    _, I2 = resolver(links2)
    assert np.abs(I1 - I2).max() < 1e-12


@pytest.mark.parametrize("nome", list(FLUXOS))
def test_estado_puro(fundos, nome):
    C, _ = fundos[nome]
    assert np.allclose(C @ C, C, atol=1e-8)


# ---------------------------------------------------------------- P2

@pytest.mark.parametrize("nome", list(FLUXOS))
def test_rede_quadrada_emerge_em_todo_fundo_de_fluxo(fundos, nome):
    _, I = fundos[nome]
    A = knn_graph(I, 4)
    assert (A & REDE).sum() / A.sum() == 1.0
    assert spearman(upper(graph_distances(A)), upper(T)) > 0.95


# ---------------------------------------------------------------- P3

def test_fluxo_muda_a_intensidade_mas_nao_a_vizinhanca(fundos):
    I0 = fundos["metal_0"][1][0, 1]
    Ipi = fundos["dirac_pi"][1][0, 1]
    assert abs(Ipi - I0) / I0 > 0.2


def test_mesma_subrede_tem_informacao_mutua_nula_a_meio_preenchimento(fundos):
    """Simetria partícula-buraco da rede bipartida: <c_i† c_j> = 0 para i, j na mesma sub-rede."""
    _, I = fundos["dirac_pi"]
    assert I[0, 2] < 1e-12 and I[0, L + 1] < 1e-12


# ---------------------------------------------------------------- P4

@pytest.mark.parametrize("forca", [0.5, np.pi])
def test_fluxo_aleatorio_nao_destroi_a_geometria(forca):
    _, I = resolver(random_flux_links(L, np.random.default_rng(2), forca))
    A = knn_graph(I, 4)
    assert (A & REDE).sum() / A.sum() > 0.95


# ---------------------------------------------------------------- P5

def test_controle_aleatorio():
    R = np.random.default_rng(0).random(T.shape)
    D = graph_distances(knn_graph(R + R.T, 4))
    assert D.max() < graph_distances(REDE).max() / 2
    assert abs(spearman(upper(D), upper(T))) < 0.1


@pytest.mark.skip(reason="Pendente: U(1) DINÂMICO em 2+1D (QED₃ compacta com férmions) exige métodos além de "
                         "diagonalização exata (DMRG em cilindro / Monte Carlo); o fundo estático é campo médio.")
def test_qed3_dinamica():
    ...
