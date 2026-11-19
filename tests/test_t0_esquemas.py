"""T0 — Dependência de esquema: o G induzido pelo Modelo Padrão é robusto ou artefato da regularização?

Predição pré-registrada (Notion): se a tese "gravidade emergente só das 3 forças + matéria do MP" for
robusta, o sinal e a ordem de grandeza de G não devem depender do regulador.
"""
from math import pi, sqrt

import pytest

from pesquisa_fisica.induced_gravity import FieldContent, inverse_G, kappa_scalar, standard_model
from pesquisa_fisica.schemes import (
    PLANCK_MASS_GEV,
    SpeciesRegulators,
    critical_vector_regulator_ratio,
    dimreg_inverse_G_SM,
    dimreg_inverse_G_massless,
    pv_inverse_G,
    pv_inverse_G_single,
    pv_weights,
)


# ---------------------------------------------------------------- Pauli-Villars: consistência

def test_regras_de_soma_PV():
    M1, M2 = 1.0, 2.3
    c1, c2 = pv_weights(M1, M2)
    assert 1 + c1 + c2 == pytest.approx(0, abs=1e-12)
    assert c1 * M1**2 + c2 * M2**2 == pytest.approx(0, abs=1e-12)


@pytest.mark.parametrize("M1,M2", [(1.0, 2.0), (3.0, 0.5), (1.0, 1.001)])
def test_PV_escalar_sempre_atrativo(M1, M2):
    """Para cada espécie, o sinal de 1/G em PV é o sinal de κ, qualquer que sejam as massas."""
    assert pv_inverse_G_single(1 / 6, M1, M2) > 0
    assert pv_inverse_G_single(-2 / 3, M1, M2) < 0


def test_PV_limite_degenerado_reproduz_corte_de_tempo_proprio():
    """Com M1 → M2 = M, PV dá κM²/(2π): a mesma forma do corte de tempo próprio com Λ = M."""
    M = 2.0
    assert pv_inverse_G_single(1 / 6, M, M * (1 + 1e-6)) == pytest.approx(
        inverse_G(FieldContent(real_scalars=1), M), rel=1e-5)


def test_PV_com_reguladores_universais_concorda_em_sinal_com_corte():
    sm = standard_model()
    assert pv_inverse_G(sm, SpeciesRegulators(1, 1, 1)) > 0


# ---------------------------------------------------------------- PV: sensibilidade do sinal

def test_razao_critica_dos_reguladores_vetoriais():
    """Matéria: 4·(1/6) + 45·(1/6) = 49/6; vetores: 12·(2/3) = 8 → razão crítica √(49/48) ≈ 1,0104."""
    r = critical_vector_regulator_ratio(standard_model())
    assert r == pytest.approx(sqrt(49 / 48), rel=1e-12)
    assert r - 1 < 0.011


def test_reguladores_vetoriais_1_porcento_mais_pesados_invertem_o_sinal():
    sm = standard_model()
    assert pv_inverse_G(sm, SpeciesRegulators(1, 1, 1.00)) > 0
    assert pv_inverse_G(sm, SpeciesRegulators(1, 1, 1.02)) < 0


def test_reguladores_vetoriais_mais_leves_reforcam_atracao():
    sm = standard_model()
    assert pv_inverse_G(sm, SpeciesRegulators(1, 1, 0.5)) > pv_inverse_G(sm, SpeciesRegulators(1, 1, 1))


# ---------------------------------------------------------------- regularização dimensional / zeta

def test_dimreg_campos_sem_massa_nao_induzem_G():
    assert dimreg_inverse_G_massless(standard_model()) == 0.0


def test_dimreg_MP_massivo_fica_33_ordens_abaixo_de_planck():
    """Só com as massas do MP, |1/G| induzido ~ m_top² ≪ M_Pl²."""
    for mu in (10.0, 173.0, 1e3):
        ratio = abs(dimreg_inverse_G_SM(mu)) / PLANCK_MASS_GEV**2
        assert ratio < 1e-30


def test_dimreg_G_depende_de_mu():
    """A parte finita corre com μ: sem contratermo não há predição de G."""
    assert dimreg_inverse_G_SM(10.0) != pytest.approx(dimreg_inverse_G_SM(1e3), rel=1e-3)


# ---------------------------------------------------------------- veredito (registrado como teste)

@pytest.mark.xfail(strict=True, reason=(
    "Veredito v0.2: o sinal de G induzido pelo MP NÃO é independente de esquema. "
    "Basta os reguladores vetoriais serem ~1% mais pesados para G < 0; em dim-reg, G não é predito."))
def test_sinal_de_G_independente_de_esquema():
    sm = standard_model()
    sinais = {
        inverse_G(sm, 1.0) > 0,
        pv_inverse_G(sm, SpeciesRegulators(1, 1, 1)) > 0,
        pv_inverse_G(sm, SpeciesRegulators(1, 1, 1.05)) > 0,
        pv_inverse_G(sm, SpeciesRegulators(1.05, 1, 1)) > 0,
    }
    assert sinais == {True}
