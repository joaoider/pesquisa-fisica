"""H2 · T1 — Z2 com matéria em REDES GRANDES: modelo de Kitaev no favo de mel (solução exata).

Kitaev (2006) = teoria de gauge Z2 (u_jk nas ligações) + matéria de Majorana, livre em cada setor de fluxo.
Exato para centenas/milhares de sítios. A parte de gauge do emaranhamento é só termo de fronteira
(Yao & Qi 2010), então a informação mútua entre células vem da matéria.

Predições pré-registradas (Notion):
  P1. O solver reproduz a energia conhecida (−0,7872 por spin, isotrópico) por duas rotas independentes.
  P2. A informação mútua entre células é invariante de gauge.
  P3. Na fase gapless B, o grafo dos 6 parceiros mais emaranhados É a rede real (precisão 1), em todos
      os tamanhos até L = 32 (1024 células, 2048 Majoranas), com o mesmo diâmetro e a mesma dimensão.
  P4. Na fase gapped A a geometria continua 2D, mas a métrica emergente é deformada pela anisotropia.
  P5. Um fundo de gauge desordenado (fluxos aleatórios) degrada, mas não destrói, a geometria.
  P6. Controle: informação "mútua" aleatória não gera geometria (diâmetro de mundo pequeno, ρ ≈ 0).
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.geometry_from_entanglement import (
    ball_growth_dimension,
    graph_distances,
    knn_graph,
    spearman,
    upper,
)
from pesquisa_fisica.emergence.kitaev_honeycomb import Honeycomb, cell_mutual_information, gauge_transform

L = 20


def resolver(L, J=(1, 1, 1), u=None):
    lat = Honeycomb(L, J=J, u=u)
    E, G = lat.ground_state_covariance()
    return lat, E, G, cell_mutual_information(G, lat.n_cells)


@pytest.fixture(scope="module")
def rede_real():
    T = Honeycomb(L).torus_cell_distances()
    A = np.isclose(T, 1.0)
    return T, A, graph_distances(A)


@pytest.fixture(scope="module")
def fases():
    rng = np.random.default_rng(3)
    return {
        "B": resolver(L, (1, 1, 1)),
        "B_anisotropica": resolver(L, (0.8, 0.6, 1)),
        "A": resolver(L, (0.3, 0.3, 1)),
        "B_fluxo_aleatorio": resolver(L, (1, 1, 1), rng.choice([-1.0, 1.0], size=3 * L * L)),
    }


# ---------------------------------------------------------------- P1: o solver está certo

def test_energia_duas_rotas_e_literatura():
    lat, E, G, _ = resolver(20)
    k = 2 * np.pi * np.arange(20) / 20
    f = np.abs(1 + np.exp(1j * k)[:, None] + np.exp(1j * k)[None, :])
    assert E / (2 * 20 * 20) == pytest.approx(-f.mean() / 2, abs=1e-10)       # rota no espaço de momentos
    assert E / (2 * 20 * 20) == pytest.approx(-0.7872, abs=5e-4)               # valor conhecido (L → ∞)
    assert np.allclose(G @ G, -np.eye(len(G)), atol=1e-8)                      # estado puro


def test_sem_vortices_tem_menor_energia_teorema_de_lieb(fases):
    assert fases["B"][1] < fases["B_fluxo_aleatorio"][1]


# ---------------------------------------------------------------- P2: invariância de gauge

def test_informacao_mutua_invariante_de_gauge():
    lat, E1, _, I1 = resolver(14)
    rng = np.random.default_rng(7)
    sitios = rng.choice(2 * lat.n_cells, size=lat.n_cells // 2, replace=False)
    _, E2, _, I2 = resolver(14, u=gauge_transform(lat, sitios))
    assert E1 == pytest.approx(E2, abs=1e-9)
    assert np.allclose(I1, I2, atol=1e-10)


# ---------------------------------------------------------------- P3: fase B reconstrói a rede exatamente

@pytest.mark.parametrize("fase", ["B", "B_anisotropica"])
def test_fase_B_grafo_emergente_e_a_rede_real(fases, rede_real, fase):
    _, _, _, I = fases[fase]
    T, A_true, D_true = rede_real
    A = knn_graph(I, 6)
    assert (A & A_true).sum() / A.sum() == 1.0
    D = graph_distances(A)
    assert D.max() == D_true.max()
    assert spearman(upper(D), upper(T)) > 0.98


@pytest.mark.parametrize("Lg", [14, 26, 32])
def test_fase_B_em_varios_tamanhos(Lg):
    lat, _, _, I = resolver(Lg)
    T = lat.torus_cell_distances()
    A_true = np.isclose(T, 1.0)
    A = knn_graph(I, 6)
    assert (A & A_true).sum() / A.sum() == 1.0
    assert ball_growth_dimension(graph_distances(A)) == pytest.approx(
        ball_growth_dimension(graph_distances(A_true)), abs=1e-9)


# ---------------------------------------------------------------- P4: fase A — 2D, métrica deformada

def test_fase_A_continua_2d_com_metrica_deformada(fases, rede_real):
    _, _, _, I = fases["A"]
    T, A_true, D_true = rede_real
    A = knn_graph(I, 6)
    assert np.all(A.sum(1) == 6)                                    # 6 vizinhos emergentes por célula
    D = graph_distances(A)
    assert np.isfinite(D).all()                                      # conectado
    assert ball_growth_dimension(D) == pytest.approx(ball_growth_dimension(D_true), abs=0.3)
    precisao = (A & A_true).sum() / A.sum()
    assert precisao == pytest.approx(2 / 3, abs=1e-9)                # 4 das 6 arestas coincidem


def test_fase_A_dimeros_fortes_encurtam_a_distancia(fases):
    """Ligações z fortes (dímeros) ligam a célula r a r ± a2: I nessa direção ≫ nas outras."""
    lat, _, _, I = fases["A"]
    ao_longo_de_a2 = I[0, lat.cell(0, 1)]
    ao_longo_de_a1 = I[0, lat.cell(1, 0)]
    assert ao_longo_de_a2 / ao_longo_de_a1 > 10


# ---------------------------------------------------------------- P5: fundo de gauge desordenado

def test_fluxo_aleatorio_degrada_mas_nao_destroi(fases, rede_real):
    _, _, _, I = fases["B_fluxo_aleatorio"]
    T, A_true, _ = rede_real
    A = knn_graph(I, 6)
    precisao = (A & A_true).sum() / A.sum()
    assert 0.6 < precisao < 1.0
    assert spearman(upper(graph_distances(A)), upper(T)) > 0.8


# ---------------------------------------------------------------- P6: controle negativo

def test_controle_aleatorio_nao_gera_geometria(rede_real):
    T, _, D_true = rede_real
    rng = np.random.default_rng(0)
    R = rng.random(T.shape)
    D = graph_distances(knn_graph(R + R.T, 6))
    assert D.max() < D_true.max() / 2                                # "mundo pequeno", sem geometria
    assert abs(spearman(upper(D), upper(T))) < 0.1
