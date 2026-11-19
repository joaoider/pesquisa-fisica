"""H2 · T1 — U(1) DINÂMICO com férmions: modelo de Schwinger na rede (QED em 1+1D), diagonalização exata.

Predições pré-registradas (Notion):
  P1. Solver correto: em J = 0 reproduz férmions livres; carga total nula.
  P2. Depois de eliminar o campo de gauge pela lei de Gauss, o Hamiltoniano de qubits é NÃO local
      (acopla qubits distantes com força que não decai com a distância).
  P3. Mesmo assim, a geometria emergente do emaranhamento é a cadeia 1D local (vizinhos recíprocos = cadeia).
  P4. Quanto mais forte o acoplamento de gauge J (confinamento), MAIS local fica o emaranhamento.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.geometry_from_entanglement import knn_graph, spearman, upper
from pesquisa_fisica.emergence.schwinger import (
    expectation_charge,
    free_fermion_energy,
    ground_state,
    zz_couplings,
)
from pesquisa_fisica.emergence.z2_gauge import reduced_density_matrix, von_neumann

N = 20                     # 20 sítios staggered = 10 células físicas; setor de carga zero: 184 756 estados
NC = N // 2
DIST = np.abs(np.subtract.outer(np.arange(NC), np.arange(NC))).astype(float)
CADEIA = np.isclose(DIST, 1.0)


def cell_mutual_information(psi, n):
    cells = [[2 * i, 2 * i + 1] for i in range(n // 2)]
    s1 = [von_neumann(reduced_density_matrix(psi, n, c)) for c in cells]
    I = np.zeros((n // 2, n // 2))
    for i in range(n // 2):
        for j in range(i + 1, n // 2):
            s12 = von_neumann(reduced_density_matrix(psi, n, cells[i] + cells[j]))
            I[i, j] = I[j, i] = max(s1[i] + s1[j] - s12, 0.0)
    return I


ACOPLAMENTOS = [0.0, 0.5, 2.0]


@pytest.fixture(scope="module")
def estados():
    out = {}
    for J in ACOPLAMENTOS:
        E, psi, _ = ground_state(N, w=1.0, m=0.5, J=J)
        out[J] = (E, psi, cell_mutual_information(psi, N))
    return out


# ---------------------------------------------------------------- P1

def test_limite_sem_gauge_reproduz_fermions_livres(estados):
    assert estados[0.0][0] == pytest.approx(free_fermion_energy(N, 1.0, 0.5), abs=1e-9)


@pytest.mark.parametrize("J", ACOPLAMENTOS)
def test_carga_total_nula(estados, J):
    assert expectation_charge(estados[J][1], N).sum() == pytest.approx(0.0, abs=1e-10)


# ---------------------------------------------------------------- P2

def test_hamiltoniano_de_qubits_e_nao_local():
    C = zz_couplings(N, J=1.0)
    n_pares = N * (N - 1) // 2
    assert np.count_nonzero(np.triu(C, 1)) >= 0.9 * n_pares              # ≥ 90% dos pares acoplados
    assert C[0, N // 2] > 5 * C[N - 3, N - 2]                            # qubits a 10 sítios acoplam 9× mais que vizinhos


# ---------------------------------------------------------------- P3

@pytest.mark.parametrize("J", ACOPLAMENTOS)
def test_geometria_emergente_e_a_cadeia_local(estados, J):
    I = estados[J][2]
    A = knn_graph(I, 2, mutual=True)
    assert (A & CADEIA).sum() == A.sum() == CADEIA.sum()                # precisão = recall = 1
    assert spearman(upper(-I), upper(DIST)) > 0.95


# ---------------------------------------------------------------- P4

def test_conclusoes_valem_na_particao_fermionica(estados):
    """Robustez: com a partição fermiônica correta (sinais de Jordan-Wigner), em vez da de qubits,
    a cadeia continua sendo recuperada e a localidade continua crescendo com J."""
    from pesquisa_fisica.emergence.su2_lattice import site_mutual_information as fermionic_mi
    from pesquisa_fisica.emergence.schwinger import sector_basis
    basis = sector_basis(N)
    r = []
    for J in ACOPLAMENTOS:
        I = fermionic_mi(basis, estados[J][1][basis], N, modes_per_site=2)
        A = knn_graph(I, 2, mutual=True)
        assert (A & CADEIA).sum() == A.sum() == CADEIA.sum()
        r.append(I[4, 5] / I[3, 6])
    assert r[0] < r[1] < r[2]


def test_acoplamento_de_gauge_mais_forte_deixa_o_emaranhamento_mais_local(estados):
    def razao(J):
        I = estados[J][2]
        return I[4, 5] / I[3, 6]                                         # vizinhos / distância 3
    r = [razao(J) for J in ACOPLAMENTOS]
    assert r[0] < r[1] < r[2]
    assert r[2] / r[0] > 100
