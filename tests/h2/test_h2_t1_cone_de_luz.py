"""H2 · T1 — Cone de luz, estrutura causal (assinatura lorentziana) e isotropia de Lorentz emergente.

Predições pré-registradas (Notion):
  P1. Em matéria local na rede, a influência de uma perturbação tem velocidade máxima finita (Lieb-Robinson):
      fronteira linear r = v t; fora do cone, amplitude exponencialmente pequena.
  P2. Estrutura causal: praticamente todo o peso fica em s² = (vt)² − r² ≥ 0 — tempo e espaço entram com
      sinais opostos (assinatura lorentziana), não como numa geometria euclidiana.
  P3. O mesmo vale numa cadeia INTERAGENTE (Heisenberg), com velocidade finita.
  P4. Isotropia de Lorentz emergente: no cone de Dirac (fluxo π, campo médio da QED₃), E = v|q| é isotrópico,
      com anisotropia ∝ q² → 0 em baixa energia. Um metal (Φ = 0) não tem isso (anisotropia O(1),
      dispersão quadrática no fundo da banda: sem velocidade-limite).
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.light_cone import (
    arrival_times,
    commutator_norms,
    dirac_anisotropy,
    fermi_velocity_anisotropy_zero_flux,
    pi_flux_dispersion,
    propagator_front,
    zero_flux_dispersion,
)

N = 401
TEMPOS = np.linspace(0, 60, 601)


@pytest.fixture(scope="module")
def frente():
    return propagator_front(N, TEMPOS)


# ---------------------------------------------------------------- P1

def test_frente_linear_com_velocidade_de_lieb_robinson(frente):
    arr = arrival_times(frente, TEMPOS, 1e-3)
    r = np.arange(len(arr))
    ok = (r >= 10) & (r <= 110) & ~np.isnan(arr)
    v = 1 / np.polyfit(r[ok], arr[ok], 1)[0]
    assert 1.9 < v < 2.3                      # velocidade máxima de grupo 2w = 2


def test_fora_do_cone_amplitude_desprezivel(frente):
    it = np.argmin(abs(TEMPOS - 20))
    assert frente[it, 60] < 1e-6              # r = 60 > v t = 40
    assert frente[it, 30] > 1e-2              # dentro do cone


# ---------------------------------------------------------------- P2

def test_estrutura_causal_quase_todo_peso_em_s2_positivo(frente):
    T, R = np.meshgrid(TEMPOS, np.arange(frente.shape[1]), indexing="ij")
    fora = R > 2 * T + 10                     # s² < 0 com margem
    assert (frente[fora] ** 2).sum() / (frente**2).sum() < 1e-4


# ---------------------------------------------------------------- P3

def test_cadeia_interagente_tem_cone_de_luz():
    ts = np.linspace(0, 3, 61)
    C = commutator_norms(8, ts)
    arr = arrival_times(C, ts, 0.05)
    assert np.all(np.diff(arr) >= 0)                      # chegada cresce com a distância
    v = 1 / np.polyfit(np.arange(2, 8), arr[2:], 1)[0]
    assert 2 < v < 20                                     # velocidade finita
    assert C[np.argmin(abs(ts - 0.3)), 7] < 1e-3          # sítio distante ainda não "sabe"


# ---------------------------------------------------------------- P4

def test_isotropia_de_lorentz_emerge_no_cone_de_dirac():
    a = [dirac_anisotropy(q) for q in (0.2, 0.1, 0.05)]
    assert a[0] / a[1] == pytest.approx(4, rel=0.05)      # ∝ q²
    assert a[1] / a[2] == pytest.approx(4, rel=0.05)
    assert a[2] < 1e-3


def test_dirac_tem_velocidade_limite_metal_nao():
    q = np.array([0.1, 0.05, 0.025])
    v_dirac = pi_flux_dispersion(np.pi / 2 + q, np.pi / 2) / q
    assert np.allclose(v_dirac, 2.0, rtol=2e-3)           # E = v|q|, v = 2t
    inclinacao_metal = (zero_flux_dispersion(q, 0) - zero_flux_dispersion(0, 0)) / q
    assert np.all(inclinacao_metal < 0.11) and inclinacao_metal[-1] < inclinacao_metal[0] / 3   # E ∝ q²


def test_metal_e_anisotropico():
    assert fermi_velocity_anisotropy_zero_flux(mu=-1.0) > 0.3
