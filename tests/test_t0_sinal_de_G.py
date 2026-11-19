"""T0 — Sinal e magnitude de G induzido pelo conteúdo do Modelo Padrão.

Predição pré-registrada (Notion, FIS-TESTE-2): 1/G > 0 e Λ ~ M_Pl.
Este arquivo documenta o RESULTADO, inclusive a fragilidade dele.
"""
from fractions import Fraction

import pytest

from pesquisa_fisica.induced_gravity import (
    FieldContent,
    cutoff_for_newton,
    standard_model,
)


def test_contagem_de_campos_do_modelo_padrao():
    sm = standard_model()
    assert (sm.real_scalars, sm.weyl_fermions, sm.gauge_vectors) == (4, 45, 12)


def test_soma_kappa_MP_acoplamento_minimo():
    """Σκ = 4(1/6) + 45(1/6) - 12(2/3) = 1/6  → em unidades de escalar: N0 + N½ - 4N1 = 4 + 45 - 48 = 1."""
    assert standard_model(xi=0).kappa_total() == Fraction(1, 6)


def test_MP_minimo_induz_gravidade_atrativa():
    assert standard_model(xi=0).kappa_total() > 0


def test_corte_necessario_e_da_ordem_de_planck():
    lam = cutoff_for_newton(standard_model(xi=0))
    assert lam == pytest.approx((12 * 3.141592653589793) ** 0.5, rel=1e-9)  # ≈ 6.14 M_Pl
    assert 1 < lam < 10


def test_xi_critico_do_higgs():
    """O sinal de G vira se o acoplamento do Higgs à curvatura passar de ξ_c = 1/24."""
    assert standard_model().critical_xi() == Fraction(1, 24)
    assert standard_model(xi=Fraction(1, 6)).kappa_total() < 0  # Higgs conforme → G < 0


def test_neutrinos_de_mao_direita_ampliam_margem():
    sm_nu = standard_model(right_handed_neutrinos=True)
    assert sm_nu.kappa_total() == Fraction(2, 3)
    assert sm_nu.critical_xi() == Fraction(1, 6)


def test_so_gauge_sem_materia_da_gravidade_repulsiva():
    """Os bósons de gauge sozinhos (as '3 forças' sem matéria) dão 1/G < 0."""
    assert FieldContent(gauge_vectors=12).kappa_total() < 0
    with pytest.raises(ValueError):
        cutoff_for_newton(FieldContent(gauge_vectors=12))
