"""Ponte H1 ↔ H2 — a entropia de área de emaranhamento é S = A/(4G) com o MESMO G induzido de H1?

Mesmo esquema de regularização (corte de tempo próprio) nas duas grandezas.

Predições pré-registradas (Notion):
  P1. O heat kernel no cone tem o termo da ponta (1/12)(1/α − α) — verificado pelo espectro exato.
  P2. Escalar mínimo: S/A pelo método das réplicas = 1/(4G_ind) de H1, exatamente (Λ²/48π).
  P3. Férmions: peso de entropia = peso em 1/G (carga central de Dirac = de escalar em 1+1D, e a
      multiplicidade espinorial transversa dá Weyl = 1 unidade de escalar).
  P4. Vetores de gauge: a contagem de estados dá +2 por vetor, mas 1/G dá −4. A ponte só fecha com o termo
      de contato de Kabat (−6 por vetor) = modos de borda do campo de gauge (Donnelly-Wall).
  P5. Modelo Padrão: sem modos de borda S_emaranhamento = 73 × A/(4G_ind); com eles, razão exatamente 1.
"""
import numpy as np
import pytest

from pesquisa_fisica.bridge import (
    entanglement_entropy_per_area_scalar,
    expected_cone_constant,
    fit_cone_constant,
    g_weight,
    standard_model_tally,
)
from pesquisa_fisica.induced_gravity import KAPPA_VECTOR, KAPPA_WEYL, FieldContent, inverse_G, kappa_scalar


# ---------------------------------------------------------------- P1

@pytest.mark.parametrize("n", [1, 2, 3, 4])
def test_heat_kernel_no_cone_pelo_espectro_exato(n):
    a = 1 / n
    f = fit_cone_constant(n)
    assert f["area"] == pytest.approx(a / 4, abs=1e-4)
    assert f["perimetro"] == pytest.approx(-a * np.sqrt(np.pi) / 4, abs=1e-4)
    assert f["constante"] == pytest.approx(expected_cone_constant(a), abs=1e-3)   # inclui (1/12)(1/α − α)


# ---------------------------------------------------------------- P2

@pytest.mark.parametrize("cutoff", [0.5, 3.7, 6.14])
def test_escalar_entropia_de_area_usa_o_G_induzido(cutoff):
    S_por_area = entanglement_entropy_per_area_scalar(cutoff)
    A_sobre_4G = inverse_G(FieldContent(real_scalars=1), cutoff) / 4
    assert S_por_area == pytest.approx(A_sobre_4G, rel=1e-8)
    assert S_por_area == pytest.approx(cutoff**2 / (48 * np.pi), rel=1e-8)


def test_pesos_de_G_coincidem_com_H1():
    k0 = 1 / 6
    assert g_weight("escalar") == pytest.approx(float(kappa_scalar(0)) / k0)
    assert g_weight("weyl") == pytest.approx(float(KAPPA_WEYL) / k0)
    assert g_weight("vetor") == pytest.approx(float(KAPPA_VECTOR) / k0)


# ---------------------------------------------------------------- P3

def test_carga_central_de_dirac_igual_a_de_escalar_em_1mais1():
    from pesquisa_fisica.emergence.entanglement import fermion_entropy, harmonic_entropy
    N = 400
    ls = np.arange(10, 200, 10)
    chord = np.log(2 * N / np.pi * np.sin(np.pi * ls / N))
    K = 2 * np.eye(N) - np.eye(N, k=1) - np.eye(N, k=-1)
    w, V = np.linalg.eigh(K)
    X, P = 0.5 * (V / np.sqrt(w)) @ V.T, 0.5 * (V * np.sqrt(w)) @ V.T
    c_escalar = 6 * np.polyfit(chord, [harmonic_entropy(X, P, list(range(l))) for l in ls], 1)[0]
    e, U = np.linalg.eigh(-(np.eye(N, k=1) + np.eye(N, k=-1)))
    C = U[:, : N // 2] @ U[:, : N // 2].T
    c_fermion = 6 * np.polyfit(chord, [fermion_entropy(C, list(range(l))) for l in ls], 1)[0]
    assert c_escalar == pytest.approx(1.0, abs=0.05)
    assert c_fermion == pytest.approx(c_escalar, abs=0.06)
    # ⇒ em 4D: Dirac = (cone 2D, c = 1) × (2 componentes espinoriais transversas) = 2 escalares; Weyl = 1
    assert g_weight("weyl") == 1.0


# ---------------------------------------------------------------- P4

def test_acoplamento_nao_minimo_xi_e_termo_de_wald():
    """No espaço plano o estado não depende de ξ (o termo ξR é nulo): a entropia por contagem de estados não muda.
    Mas 1/G depende de ξ: o escalar conforme (ξ = 1/6) tem entropia de emaranhamento e NÃO induz G."""
    assert g_weight("escalar", xi=1 / 6) == pytest.approx(0.0)


def test_vetores_contagem_de_estados_contra_G_induzido():
    assert g_weight("vetor") == -4.0           # contribui contra a gravidade (H1)
    # contagem de estados: 2 polarizações físicas ⇒ +2; diferença −6 = termo de contato de Kabat


# ---------------------------------------------------------------- P5

@pytest.mark.xfail(strict=True, reason="Sem modos de borda do gauge, a entropia de emaranhamento do MP é 73× A/(4G_ind): "
                                       "a ponte H1↔H2 NÃO fecha só com graus de liberdade de bulk.")
def test_ponte_fecha_sem_modos_de_borda():
    assert standard_model_tally()["razao"] == pytest.approx(1.0)


def test_mp_sem_modos_de_borda_razao_73():
    t = standard_model_tally()
    assert t["S_emaranhamento"] == 73 and t["A_sobre_4G"] == 1 and t["razao"] == 73


@pytest.mark.parametrize("xi,nu_R", [(0.0, False), (0.0, True), (0.02, False)])
def test_ponte_fecha_com_modos_de_borda(xi, nu_R):
    t = standard_model_tally(xi=xi, edge_modes=True, right_handed_neutrinos=nu_R)
    assert t["razao"] == pytest.approx(1.0)


@pytest.mark.skip(reason="Pendente: calcular o termo de contato (−6 por vetor) diretamente, via entropia de Maxwell na rede "
                         "com centro elétrico / espaço de Hilbert estendido (Casini-Huerta; Donnelly-Wall), em vez de "
                         "tomá-lo da literatura.")
def test_termo_de_contato_de_kabat_calculado_na_rede():
    ...
