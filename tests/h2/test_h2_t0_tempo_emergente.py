"""H2 · T0 — Tempo emergente (Page-Wootters).

Predição: um universo num estado global estacionário (H_total|Ψ> = 0, "sem tempo") contém, no
estado do sistema condicionado à leitura de um relógio interno, a evolução de Schrödinger exata.
E o tempo só emerge se relógio e sistema estiverem EMARANHADOS.
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.page_wootters import (
    build,
    clock_system_entanglement,
    conditional_state,
    fidelity,
    schrodinger,
)

LEVELS = np.array([0, 1, 3, 4])          # espectro do sistema em unidades de ω
PSI0 = np.array([0.5, 0.5j, -0.5, 0.5])
N_CLOCK = 32


@pytest.fixture(scope="module")
def universo():
    return build(LEVELS, PSI0, N_CLOCK)


def test_estado_global_e_estacionario(universo):
    Psi, H_tot, _ = universo
    assert np.allclose(H_tot * Psi.ravel(), 0)


@pytest.mark.parametrize("t", [0.0, 0.3, 1.7, 5.0, 12.4])
def test_estado_condicional_obedece_schrodinger(universo, t):
    Psi, _, E_c = universo
    assert fidelity(conditional_state(Psi, t, E_c), schrodinger(PSI0, LEVELS, t)) == pytest.approx(1.0, abs=1e-12)


def test_ha_evolucao_nao_trivial(universo):
    Psi, _, E_c = universo
    assert fidelity(conditional_state(Psi, 0.0, E_c), conditional_state(Psi, 1.0, E_c)) < 0.9


def test_sem_emaranhamento_nao_ha_tempo():
    """Sistema num autoestado de energia: Ψ é produto relógio⊗sistema e o estado condicional não muda."""
    psi = np.array([0, 1, 0, 0], dtype=complex)
    Psi, _, E_c = build(LEVELS, psi, N_CLOCK)
    assert clock_system_entanglement(Psi) == pytest.approx(0.0, abs=1e-12)
    assert fidelity(conditional_state(Psi, 0.0, E_c), conditional_state(Psi, 3.0, E_c)) == pytest.approx(1.0)


def test_emaranhamento_relogio_sistema_positivo_quando_ha_dinamica(universo):
    Psi, _, _ = universo
    assert clock_system_entanglement(Psi) == pytest.approx(np.log(4), abs=1e-12)  # 4 níveis igualmente pesados
