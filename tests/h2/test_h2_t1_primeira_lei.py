"""H2 · T1 — Rumo a Einstein: Hamiltoniano de emaranhamento geométrico e primeira lei do emaranhamento.

Ingredientes das derivações de Jacobson (1995, 2015) e Faulkner et al. (2014), testados em matéria de rede
(cadeia crítica de férmions livres: CFT c = 1, v = 2):

Predições pré-registradas (Notion):
  P1. O Hamiltoniano de emaranhamento do vácuo num intervalo é LOCAL e segue o perfil de Bisognano-Wichmann
      β(x) = 2π x(ℓ−x)/(ℓ v) ("temperatura de Unruh" local): R² > 0,99, coeficiente a 10% do contínuo.
  P2. Primeira lei: δS = δ⟨K⟩ em primeira ordem, para perturbações da "métrica" local (salto modulado),
      tanto com o K exato quanto com o K geométrico (BW): razão a 2% de 1, independente de ε.
  P3. Sanidade: reescalar o salto uniformemente (só muda a unidade de tempo) não muda S nem K.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.modular import (
    bw_profile,
    entanglement_hamiltonian,
    first_law_bonds,
    infinite_chain_correlation,
)


# ---------------------------------------------------------------- P1

@pytest.mark.parametrize("ell", [12, 20, 30])
def test_hamiltoniano_de_emaranhamento_segue_bisognano_wichmann(ell):
    H = entanglement_hamiltonian(infinite_chain_correlation(ell))
    nn = np.array([H[x, x + 1] for x in range(ell - 1)])
    beta = bw_profile(ell)
    a = np.dot(nn, beta) / np.dot(beta, beta)
    r2 = 1 - np.sum((nn - a * beta) ** 2) / np.sum((nn - nn.mean()) ** 2)
    assert r2 > 0.99
    assert a == pytest.approx(-1.0, abs=0.1)                 # termo de ligação −(c†c + h.c.) com peso β


@pytest.mark.parametrize("ell", [20, 30])
def test_hamiltoniano_de_emaranhamento_e_local(ell):
    H = entanglement_hamiltonian(infinite_chain_correlation(ell))
    assert np.abs(np.diag(H)).max() < 1e-10                  # sem termo on-site (simetria partícula-buraco)
    assert max(abs(H[x, x + 2]) for x in range(ell - 2)) < 1e-10     # alcance par nulo
    alcance3 = max(abs(H[x, x + 3]) for x in range(ell - 3))
    assert alcance3 < 0.05 * max(abs(H[x, x + 1]) for x in range(ell - 1))   # alcance 3: correção pequena


# ---------------------------------------------------------------- P2

N = 2000
ELL = 40


@pytest.fixture(scope="module")
def K_exato():
    return entanglement_hamiltonian(infinite_chain_correlation(ELL))


def modulacoes():
    xb = np.arange(N - 1) + 0.5
    return {
        "onda_longa": np.cos(2 * np.pi * (xb - N / 2) / (4 * ELL)),
        "gauss_central": np.exp(-((xb - N / 2) / (ELL / 4)) ** 2),
        "gauss_borda": np.exp(-((xb - N / 2 - ELL / 2) / (ELL / 6)) ** 2),
    }


@pytest.mark.parametrize("nome", ["onda_longa", "gauss_central", "gauss_borda"])
def test_primeira_lei_do_emaranhamento(K_exato, nome):
    dS, dK_bw, dK_ex = first_law_bonds(N, ELL, modulacoes()[nome], 1e-4, K_exato)
    assert dS / dK_ex == pytest.approx(1.0, abs=0.02)
    assert dS / dK_bw == pytest.approx(1.0, abs=0.02)


def test_primeira_lei_e_linear_em_epsilon(K_exato):
    f = modulacoes()["gauss_central"]
    r = [first_law_bonds(N, ELL, f, eps, K_exato) for eps in (1e-3, 1e-4)]
    assert r[0][0] / r[1][0] == pytest.approx(10, rel=0.01)                 # δS ∝ ε
    assert r[0][0] / r[0][1] == pytest.approx(r[1][0] / r[1][1], abs=0.01)  # razão independe de ε


# ---------------------------------------------------------------- P3

def test_reescala_uniforme_nao_muda_nada(K_exato):
    dS, dK_bw, dK_ex = first_law_bonds(N, ELL, np.ones(N - 1), 1e-3, K_exato)
    assert abs(dS) < 1e-10 and abs(dK_bw) < 1e-10 and abs(dK_ex) < 1e-10


@pytest.mark.skip(reason="Pendente: passo holográfico — impor δS = δ⟨K⟩ para TODAS as regiões e reconstruir a "
                         "geometria de um bulk com dimensão extra (Faulkner et al. 2014: ⇒ Einstein linearizado). "
                         "Exige um modelo de brinquedo holográfico (redes tensoriais/HaPPY ou MERA).")
def test_einstein_linearizado_holografico():
    ...
