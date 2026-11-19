"""H2 · T1 — Espaço a partir do emaranhamento (Cao, Carroll & Michalakis 2017).

Só com a informação mútua entre graus de liberdade (sem posições), reconstruir:
  (a) distâncias que respeitam a ordem das distâncias verdadeiras;
  (b) a dimensão do espaço (1, 2 ou 3);
e mostrar que um estado SEM estrutura local não produz geometria (controle negativo).

LIMITAÇÃO (registrada em teste skip abaixo): o Hamiltoniano de brinquedo já é local numa rede.
O teste mostra que o ESTADO codifica a geometria, não de onde vem a localidade.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.entanglement import harmonic_correlations, harmonic_mutual_information_matrix
from pesquisa_fisica.emergence.geometry_from_entanglement import (
    emergent_distances,
    spearman,
    true_distances,
    upper,
)


def correlation_dimension(D: np.ndarray) -> float:
    """Inclinação de ln N(r) × ln r (N = fração de pares a distância ≤ r)."""
    d = upper(D)
    rs = np.quantile(d, np.linspace(0.05, 0.4, 20))
    cnt = np.array([(d <= r).mean() for r in rs])
    return float(np.polyfit(np.log(rs), np.log(cnt), 1)[0])


CASOS = {1: (80,), 2: (14, 14), 3: (7, 7, 7)}


@pytest.fixture(scope="module")
def distancias_emergentes():
    out = {}
    for dim, shape in CASOS.items():
        X, P = harmonic_correlations(shape, mass=0.2)
        out[dim] = (emergent_distances(harmonic_mutual_information_matrix(X, P)), true_distances(shape))
    return out


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_distancias_emergentes_preservam_a_ordem(distancias_emergentes, dim):
    D, T = distancias_emergentes[dim]
    assert spearman(upper(D), upper(T)) > 0.98


@pytest.mark.parametrize("dim", [1, 2, 3])
def test_dimensao_emergente(distancias_emergentes, dim):
    D, _ = distancias_emergentes[dim]
    assert correlation_dimension(D) == pytest.approx(dim, abs=0.35)


def test_dimensoes_em_ordem_crescente(distancias_emergentes):
    dims = [correlation_dimension(distancias_emergentes[d][0]) for d in (1, 2, 3)]
    assert dims[0] < dims[1] < dims[2]


def test_controle_negativo_grafo_aleatorio_nao_gera_espaco_de_baixa_dimensao():
    """Acoplamentos não locais aleatórios: o emaranhamento existe, mas não define um espaço 2D."""
    rng = np.random.default_rng(0)
    n = 144
    A = np.triu(rng.random((n, n)) < 0.05, 1)
    A = (A | A.T).astype(float)
    K = 0.04 * np.eye(n) + np.diag(A.sum(1)) - A
    w, V = np.linalg.eigh(K)
    X, P = 0.5 * (V / np.sqrt(w)) @ V.T, 0.5 * (V * np.sqrt(w)) @ V.T
    D = emergent_distances(harmonic_mutual_information_matrix(X, P))
    assert correlation_dimension(D) > 4
    assert abs(spearman(upper(D), upper(true_distances((12, 12))))) < 0.1


def test_estado_produto_nao_define_geometria():
    with pytest.raises(ValueError):
        emergent_distances(np.zeros((5, 5)))


@pytest.mark.skip(reason="Coberto em test_h2_t0_localidade_do_espectro.py (v0.9): caso invariante por "
                         "translação passou; caso desordenado inconclusivo.")
def test_localidade_emergente_do_espectro():
    ...
