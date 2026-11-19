"""H2 · T1 — Geometria a partir de uma teoria de GAUGE na rede: Z2 com matéria (Fradkin-Shenker).

Primeiro passo de "matéria livre" → "forças". Diagonalização exata, rede 3×3 periódica (18 qubits
de ligação, invariantes de gauge no calibre unitário).

Predições pré-registradas (Notion):
  P1. Código tórico puro (h = λ = 0: só a estrutura de gauge, sem flutuações): informação mútua
      entre ligações é EXATAMENTE zero → nenhuma geometria local emerge.
  P2. Com campo elétrico (h) e/ou matéria (λ), os pares mais emaranhados são os vizinhos reais
      (precisão 1) em todas as fases: desconfinada, confinada e de Higgs.
  P3. Dualidade elétrico-magnética: (h, 0) e (0, h) dão a mesma energia e o mesmo espectro de I.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.geometry_from_entanglement import emergent_distances, upper
from pesquisa_fisica.emergence.z2_gauge import (
    SquareTorus,
    expectation_product_x,
    ground_state,
    mutual_information_matrix,
)

LAT = SquareTorus(3)
N = LAT.n_links
DIST = upper(LAT.torus_distances())
NEAR = DIST <= 1.0 + 1e-9           # vizinhos: ligações perpendiculares (0,71) e paralelas (1,0)

FASES = {
    "codigo_torico": (0.0, 0.0),
    "desconfinada_eletrico": (0.3, 0.0),
    "desconfinada_materia": (0.0, 0.3),
    "mista": (0.2, 0.2),
    "perto_da_transicao": (0.6, 0.6),
    "confinada": (1.5, 0.0),
    "higgs": (0.0, 1.5),
}


@pytest.fixture(scope="module")
def estados():
    out = {}
    for nome, (h, lam) in FASES.items():
        E, psi = ground_state(LAT, h=h, lam=lam)
        out[nome] = (E, psi, mutual_information_matrix(psi, N))
    return out


# ---------------------------------------------------------------- P1: forças puras não geram geometria

def test_codigo_torico_energia_e_estabilizadores(estados):
    E, psi, _ = estados["codigo_torico"]
    assert E == pytest.approx(-2 * LAT.L**2)
    assert all(expectation_product_x(psi, N, s) == pytest.approx(1.0) for s in LAT.stars)


def test_codigo_torico_informacao_mutua_nula(estados):
    _, _, I = estados["codigo_torico"]
    assert I.max() < 1e-12


def test_codigo_torico_nao_define_geometria(estados):
    with pytest.raises(ValueError):
        emergent_distances(estados["codigo_torico"][2])


# ---------------------------------------------------------------- P2: com flutuações, a vizinhança emerge

@pytest.mark.parametrize("fase", ["desconfinada_eletrico", "desconfinada_materia", "mista", "confinada", "higgs"])
def test_pares_mais_emaranhados_sao_os_vizinhos(estados, fase):
    u = upper(estados[fase][2])
    top = np.argsort(u)[::-1][: NEAR.sum()]
    assert NEAR[top].mean() == 1.0


@pytest.mark.parametrize("fase", ["desconfinada_eletrico", "mista", "confinada", "higgs"])
def test_vizinhos_muito_mais_emaranhados_que_distantes(estados, fase):
    u = upper(estados[fase][2])
    assert u[NEAR].mean() / u[~NEAR].mean() > 10


def test_geometria_borra_perto_da_transicao(estados):
    """Perto da transição o comprimento de correlação cresce e o contraste perto/longe cai."""
    def ratio(f):
        u = upper(estados[f][2])
        return u[NEAR].mean() / u[~NEAR].mean()
    assert ratio("perto_da_transicao") < ratio("mista") / 5


# ---------------------------------------------------------------- P3: dualidade elétrico-magnética

def test_dualidade_h_lambda(estados):
    E1, _, I1 = estados["desconfinada_eletrico"]
    E2, _, I2 = estados["desconfinada_materia"]
    assert E1 == pytest.approx(E2, abs=1e-8)
    assert np.allclose(np.sort(upper(I1)), np.sort(upper(I2)), atol=1e-8)


# ---------------------------------------------------------------- pendentes

@pytest.mark.skip(reason="Pendente: U(1) (eletromagnetismo) e SU(2)/SU(3) na rede, com matéria fermiônica.")
def test_geometria_de_u1_e_su_n():
    ...
