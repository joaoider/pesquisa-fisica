"""T1 — Teste de calibração: coeficientes de heat kernel e termo de Einstein-Hilbert induzido.

Estratégia: duas rotas independentes precisam concordar.
  (a) fórmula simbólica de Seeley-DeWitt aplicada à geometria;
  (b) espectro EXATO do Laplaciano na esfera, somado numericamente e ajustado em s → 0.
Se (a) e (b) concordam, o coeficiente de R (que vira 1/G) está calibrado.
"""
import numpy as np
import pytest
import sympy as sp

from pesquisa_fisica.geometry import schwarzschild, sphere
from pesquisa_fisica.heat_kernel import fit_heat_coefficients, scalar_a1, scalar_a2


# ------------------------------------------------ rota (a): simbólica

def test_a1_esfera_unitaria_2d():
    assert sp.simplify(scalar_a1(sphere(2, radius=1)) - sp.Rational(1, 3)) == 0


@pytest.mark.parametrize("xi", [0, sp.Rational(1, 6), sp.Rational(1, 2)])
def test_a1_depende_de_xi(xi):
    m = sphere(2, radius=1)
    assert sp.simplify(scalar_a1(m, xi) - (sp.Rational(1, 6) - xi) * 2) == 0


def test_a2_esfera_unitaria_2d():
    assert sp.simplify(scalar_a2(sphere(2, radius=1)) - sp.Rational(1, 15)) == 0


def test_a1_nulo_e_a2_kretschmann_em_schwarzschild():
    m = schwarzschild()
    M, r = sp.Symbol("M", positive=True), sp.Symbol("r", real=True)
    assert scalar_a1(m) == 0
    a2 = scalar_a2(m, include_total_derivatives=False)
    assert sp.simplify(a2 - 48 * M**2 / r**6 / 180) == 0


# ------------------------------------------------ rota (b): espectro exato

@pytest.mark.parametrize("xi", [0.0, 1 / 6, 0.3])
def test_espectro_S2_confirma_a1(xi):
    b = fit_heat_coefficients(dim=2, xi=xi)
    assert b[0] == pytest.approx(1.0, abs=1e-6)
    assert b[1] == pytest.approx((1 / 6 - xi) * 2, abs=1e-5)


def test_espectro_S2_confirma_a2():
    b = fit_heat_coefficients(dim=2, xi=0.0)
    assert b[2] == pytest.approx(1 / 15, abs=1e-3)


def test_espectro_S4_confirma_a1_em_quatro_dimensoes():
    """Em S^4 unitária, R = 12: a_1 = R/6 = 2. É o coeficiente que vira 1/G em 4D."""
    b = fit_heat_coefficients(dim=4, xi=0.0, s_min=0.004, s_max=0.06)
    assert b[0] == pytest.approx(1.0, abs=1e-5)
    assert b[1] == pytest.approx(2.0, abs=1e-3)


def test_espectro_S4_confirma_a2():
    """a_2 em S^4: (Riem² - Ric²)/180 + R²/72 = (24 - 36)/180 + 144/72 = 29/15."""
    b = fit_heat_coefficients(dim=4, xi=0.0, s_min=0.004, s_max=0.06)
    assert b[2] == pytest.approx(29 / 15, rel=2e-2)


# ------------------------------------------------ termo de Einstein-Hilbert induzido

def test_G_induzido_por_um_escalar_minimo():
    """Γ ⊃ -(Λ²/32π²)·(1/6)∫√g R  ⇒  1/G = Λ²/(12π)."""
    from pesquisa_fisica.induced_gravity import FieldContent, inverse_G

    lam = 3.7
    assert inverse_G(FieldContent(real_scalars=1), lam) == pytest.approx(lam**2 / (12 * np.pi))


def test_escalar_conforme_nao_induz_G():
    from pesquisa_fisica.induced_gravity import FieldContent

    assert FieldContent(real_scalars=1, xi=sp.Rational(1, 6)).kappa_total() == 0
