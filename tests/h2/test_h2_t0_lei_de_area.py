"""H2 · T0 — O emaranhamento da matéria tem a estrutura que uma geometria exige?

Pré-condição para "espaço a partir do emaranhamento": a entropia de uma região deve escalar com a
FRONTEIRA (lei de área), não com o volume, como a entropia de Bekenstein-Hawking. Em 1D crítico,
a assinatura é S = (c/3) ln L (Calabrese-Cardy), com c = 1 para férmions livres.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.entanglement import (
    fermion_chain_correlations,
    fermion_entropy,
    harmonic_correlations,
    harmonic_entropy,
)


def test_carga_central_de_fermions_livres_e_1():
    N = 400
    C = fermion_chain_correlations(N)
    Ls = np.arange(8, 120, 8)
    S = [fermion_entropy(C, list(range(L))) for L in Ls]
    chord = np.log(N / np.pi * np.sin(np.pi * Ls / N))
    c = 3 * np.polyfit(chord, S, 1)[0]
    assert c == pytest.approx(1.0, abs=0.01)


@pytest.fixture(scope="module")
def rede_2d():
    X, P = harmonic_correlations((30, 30), mass=0.5, periodic=True)
    return X, P, np.arange(900).reshape(30, 30)


def test_lei_de_area_2d_entropia_linear_no_perimetro(rede_2d):
    X, P, idx = rede_2d
    ls = np.array([2, 3, 4, 5, 6, 8, 10])
    S = np.array([harmonic_entropy(X, P, idx[:l, :l].ravel()) for l in ls])
    r = np.corrcoef(4 * ls, S)[0, 1]
    assert r > 0.9999


def test_lei_de_volume_excluida(rede_2d):
    """S/ℓ² cai como 1/ℓ: a entropia não cresce com a área (volume em 2D) da região."""
    X, P, idx = rede_2d
    s_small = harmonic_entropy(X, P, idx[:3, :3].ravel()) / 9
    s_large = harmonic_entropy(X, P, idx[:10, :10].ravel()) / 100
    assert s_large < 0.4 * s_small


def test_estado_produto_nao_tem_emaranhamento():
    """Controle: massa enorme → sítios desacoplados → S ≈ 0. Sem emaranhamento, nada a reconstruir."""
    X, P = harmonic_correlations((8, 8), mass=1e4)
    assert harmonic_entropy(X, P, list(range(16))) < 1e-6
