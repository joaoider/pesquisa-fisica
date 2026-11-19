"""H2 · T1 — SU(3) na rede em 1+1D com férmions (QCD em 1+1D, 3 cores), exato.

Versão rápida (N = 6 sítios × 3 cores = 18 modos; 48 620 estados), parte da suíte normal.
A versão grande (N = 8, 2,7 milhões de estados) roda no GitHub Actions: scripts/su3_n8.py,
workflow "Cálculo pesado" (resultado em results/su3_n8.json).

Predições pré-registradas (Notion), as mesmas do SU(2):
  P1. Solver correto: código geral SU(N) reproduz o SU(2) já testado; J = 0 = 3 cópias de férmions livres;
      estado fundamental é singleto de cor.
  P2. Geometria emergente do emaranhamento = cadeia local, para todo J.
  P3. Confinamento de cor mais forte ⇒ emaranhamento mais local (razão vizinhos / distância 3 cresce).
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.geometry_from_entanglement import knn_graph, spearman, upper
from pesquisa_fisica.emergence.su2_lattice import SU2Chain
from pesquisa_fisica.emergence.sun_lattice import SUNChain, site_mutual_information

N = 6
DIST = np.abs(np.subtract.outer(np.arange(N), np.arange(N))).astype(float)
CADEIA = np.isclose(DIST, 1.0)
ACOPLAMENTOS = [0.0, 0.5, 2.0, 5.0]


@pytest.fixture(scope="module")
def estados():
    out = {}
    for J in ACOPLAMENTOS:
        ch = SUNChain(N, Nc=3, w=1.0, m=0.5, J=J)
        E, v = ch.ground_state()
        out[J] = (ch, E, v, site_mutual_information(ch, v))
    return out


# ---------------------------------------------------------------- P1

@pytest.mark.parametrize("J", [0.0, 1.3])
def test_codigo_geral_reproduz_su2(J):
    E_geral, _ = SUNChain(6, Nc=2, w=1.0, m=0.5, J=J).ground_state()
    E_su2, _ = SU2Chain(6, w=1.0, m=0.5, J=J).ground_state()
    assert E_geral == pytest.approx(E_su2, abs=1e-9)


def test_limite_sem_gauge_reproduz_tres_copias_de_fermions_livres(estados):
    ch, E, _, _ = estados[0.0]
    assert E == pytest.approx(ch.free_energy(), abs=1e-9)


@pytest.mark.parametrize("J", ACOPLAMENTOS)
def test_estado_fundamental_e_singleto_de_cor(estados, J):
    ch, _, v, _ = estados[J]
    assert ch.total_color_casimir(v) == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------- P2

@pytest.mark.parametrize("J", ACOPLAMENTOS)
def test_geometria_emergente_e_a_cadeia_local(estados, J):
    I = estados[J][3]
    A = knn_graph(I, 2, mutual=True)
    assert (A & CADEIA).sum() == A.sum() == CADEIA.sum()


@pytest.mark.parametrize("J", [0.5, 2.0, 5.0])
def test_com_forca_ligada_emaranhamento_ordena_as_distancias(estados, J):
    assert spearman(upper(-estados[J][3]), upper(DIST)) > 0.9


def test_sem_forca_a_ordem_por_distancia_se_inverte(estados):
    """Achado (v0.8): com J = 0 (férmions livres de 3 cores) a informação mútua a distância 3 SUPERA a de
    distância 2 (as sub-redes se correlacionam entre si): a matéria sozinha não define uma métrica
    monotônica. Qualquer J > 0 restaura a ordem. Em N = 4 o efeito chega a ligar as pontas da cadeia."""
    def I_por_distancia(I, d):
        return np.mean([I[i, i + d] for i in range(N - d)])
    I0 = estados[0.0][3]
    assert I_por_distancia(I0, 3) > I_por_distancia(I0, 2)
    assert spearman(upper(-I0), upper(DIST)) < 0.8
    for J in (0.5, 2.0, 5.0):
        I = estados[J][3]
        assert I_por_distancia(I, 2) > I_por_distancia(I, 3)


# ---------------------------------------------------------------- P3

def test_confinamento_de_cor_deixa_o_emaranhamento_mais_local(estados):
    r = [estados[J][3][2, 3] / estados[J][3][1, 4] for J in ACOPLAMENTOS]   # vizinhos / distância 3
    assert all(a < b for a, b in zip(r, r[1:]))
    assert r[-1] / r[0] > 20
