"""H2 · T1 — Teste discriminante: gap de gauge × gap sem gauge (mesmo gap).

Predições pré-registradas (Notion, 19/11/2026):
  H0 (literatura): com o mesmo gap, o perfil de informação mútua I(d) é o mesmo, venha o gap da força de gauge
      ou de uma massa sem gauge ("o gap é tudo").
  H1 (algo específico da força): perfis diferentes.
Resultado: H0 REJEITADA em U(1) e SU(2) em 1+1D (redes pequenas). No mesmo gap, a força confinante dá
emaranhamento MAIS LOCAL; o excesso cresce com o acoplamento. Mecanismo provável: bosonização — na teoria
confinada as excitações com gap são mésons; na livre, o portador mais leve é um férmion com metade do gap.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.discriminant import (
    decay_rate,
    match_mass_to_gap,
    schwinger_gap,
    schwinger_gap_and_profile,
    su2_gap_and_profile,
)

N = 16


@pytest.fixture(scope="module")
def schwinger_pareado():
    out = {}
    for J in (0.5, 1.0, 2.0):
        g, pg = schwinger_gap_and_profile(N, 0.0, J)              # gap só da força (massa de quark zero)
        m = match_mass_to_gap(lambda mm: schwinger_gap(N, mm, 0.0), g, iters=24)
        gm, pm = schwinger_gap_and_profile(N, m, 0.0)            # mesmo gap, sem gauge
        out[J] = dict(gap=g, gap_livre=gm, m=m, gauge=pg, livre=pm)
    return out


def test_gaps_pareados(schwinger_pareado):
    for r in schwinger_pareado.values():
        assert r["gap"] == pytest.approx(r["gap_livre"], abs=1e-5)


@pytest.mark.xfail(strict=True, reason="H0 rejeitada: no mesmo gap, os perfis de informação mútua NÃO coincidem.")
def test_hipotese_nula_mesmo_gap_mesmo_perfil(schwinger_pareado):
    r = schwinger_pareado[2.0]
    assert np.allclose(r["gauge"][:5], r["livre"][:5], rtol=0.2)


def test_mesmo_gap_forca_da_emaranhamento_mais_local(schwinger_pareado):
    for J, r in schwinger_pareado.items():
        assert r["gauge"][4] < r["livre"][4]                      # I(d = 5) menor com gauge
        assert decay_rate(r["gauge"], 2, 4) > decay_rate(r["livre"], 2, 4)


def test_excesso_de_localidade_cresce_com_o_acoplamento(schwinger_pareado):
    razao = [schwinger_pareado[J]["gauge"][4] / schwinger_pareado[J]["livre"][4] for J in (0.5, 1.0, 2.0)]
    assert razao[0] > razao[1] > razao[2]
    assert razao[2] < 0.2


def test_controle_mesmo_comprimento_de_correlacao_ainda_difere(schwinger_pareado):
    """Igualando a taxa de decaimento (ξ) em vez do gap, os perfis ainda diferem: a força dá MAIS informação
    mútua a curta distância. Ou seja, não é 'só ξ'."""
    alvo = decay_rate(schwinger_pareado[2.0]["gauge"], 2, 4)
    lo, hi = 0.0, 8.0
    for _ in range(22):
        mid = (lo + hi) / 2
        if decay_rate(schwinger_gap_and_profile(N, mid, 0.0)[1], 2, 4) < alvo:
            lo = mid
        else:
            hi = mid
    _, p_xi = schwinger_gap_and_profile(N, (lo + hi) / 2, 0.0)
    assert decay_rate(p_xi, 2, 4) == pytest.approx(alvo, rel=1e-3)
    assert schwinger_pareado[2.0]["gauge"][0] > 1.2 * p_xi[0]


# ---------------------------------------------------------------- SU(2)

def test_su2_mesmo_gap_tambem_mais_local():
    g, pg, S2 = su2_gap_and_profile(8, 0.5, 2.0)
    assert S2 == pytest.approx(0.0, abs=1e-8)
    m = match_mass_to_gap(lambda mm: su2_gap_and_profile(8, mm, 0.0)[0], g, iters=22)
    gm, pm, _ = su2_gap_and_profile(8, m, 0.0)
    assert gm == pytest.approx(g, abs=1e-4)
    assert pg[4] / pm[4] < 0.5


def test_su2_confinamento_sem_massa_de_quark_nao_gera_gap_crescente():
    """Confinamento ≠ gap: em SU(2) com quarks sem massa, o gap DIMINUI com o acoplamento (méson leve, à la
    't Hooft), ao contrário do U(1)."""
    gaps = [su2_gap_and_profile(8, 0.0, J)[0] for J in (0.5, 1.0, 2.0)]
    assert gaps[0] > gaps[1] > gaps[2]
