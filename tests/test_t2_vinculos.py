"""T2 — Vínculos experimentais.

Testes marcados com xfail registram problemas conhecidos e ABERTOS: falham hoje, e isso é
informação. Quando um mecanismo novo for proposto, o xfail deve virar teste normal.
"""
import pytest

from pesquisa_fisica.induced_gravity import cutoff_for_newton, induced_vacuum_energy_ratio, standard_model


@pytest.mark.xfail(strict=True, reason="Problema da constante cosmológica: a_0 induz ρ_vac ~ Λ⁴, ~10^120 acima do observado.")
def test_constante_cosmologica_induzida_compativel():
    sm = standard_model()
    lam = cutoff_for_newton(sm)
    n_dof = 4 + 2 * 12 - 2 * 45  # graus de liberdade bosônicos - fermiônicos (Weyl = 2)
    assert induced_vacuum_energy_ratio(lam, n_dof) < 10


def test_tamanho_do_problema_da_constante_cosmologica():
    sm = standard_model()
    lam = cutoff_for_newton(sm)
    ratio = induced_vacuum_energy_ratio(lam, 4 + 2 * 12 - 2 * 45)
    assert ratio > 1e120


@pytest.mark.skip(reason="Pendente: precisa da ação efetiva quadrática em h_ab (velocidade de c_g, GW170817).")
def test_velocidade_das_ondas_gravitacionais():
    ...


@pytest.mark.skip(reason="Pendente: universalidade do acoplamento e setor escuro (princípio de equivalência).")
def test_universalidade_do_acoplamento():
    ...
