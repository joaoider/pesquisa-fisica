"""T1 — Verificação independente dos pesos κ por spin (FIS-TESTE: conferir pesos contra a literatura).

Os pesos usados em induced_gravity.py vêm da fórmula de Seeley-DeWitt (Vassilevich 2003).
Aqui eles são recalculados SEM a fórmula, só com os espectros exatos de Dirac e de 1-formas na S⁴
unitária (R = 12). Se os dois caminhos concordam, os pesos estão certos (no esquema de um laço).
"""
import pytest

from pesquisa_fisica.heat_kernel import fit_spin_coefficients
from pesquisa_fisica.induced_gravity import KAPPA_VECTOR, KAPPA_WEYL

R_S4 = 12.0


def test_dirac_numero_de_componentes():
    b = fit_spin_coefficients("dirac")
    assert b[0] == pytest.approx(4.0, abs=1e-4)


def test_dirac_a1_igual_menos_R_sobre_3():
    """tr a_1(D̸²) = 4(1/6 - 1/4) R = -R/3 = -4 na S⁴ unitária."""
    b = fit_spin_coefficients("dirac")
    assert b[1] == pytest.approx(-R_S4 / 3, abs=2e-3)


def test_peso_weyl_a_partir_do_espectro():
    """Férmion: Γ = +½∫(ds/s)Tr e^{-sD̸²} (sinal oposto ao bóson) → κ_Dirac = -tr a_1/R; Weyl = metade."""
    b = fit_spin_coefficients("dirac")
    kappa_weyl_numerico = -b[1] / R_S4 / 2
    assert kappa_weyl_numerico == pytest.approx(float(KAPPA_WEYL), abs=1e-3)


def test_1formas_numero_de_componentes_e_a1():
    """Δ₁ = -∇² + Ric: E = -Ric, tr a_1 = (1/6)(-6R + 4R) = -R/3."""
    b = fit_spin_coefficients("hodge1")
    assert b[0] == pytest.approx(4.0, abs=1e-4)
    assert b[1] == pytest.approx(-R_S4 / 3, abs=2e-3)


def test_peso_vetor_de_gauge_com_fantasmas():
    b = fit_spin_coefficients("gauge_vector")
    assert b[0] == pytest.approx(2.0, abs=1e-4)  # 2 polarizações físicas
    assert b[1] / R_S4 == pytest.approx(float(KAPPA_VECTOR), abs=1e-3)


def test_combinacao_N0_mais_Nmeio_menos_4N1():
    """Em unidades de escalar (κ/κ₀): escalar 1, Weyl 1, vetor -4."""
    k0 = 1 / 6
    assert float(KAPPA_WEYL) / k0 == pytest.approx(1)
    assert float(KAPPA_VECTOR) / k0 == pytest.approx(-4)
