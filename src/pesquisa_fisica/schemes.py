"""Dependência de esquema do G induzido.

Todos os esquemas partem do mesmo termo a_1 da ação efetiva de um laço:
    Γ ⊃ -(1/32π²) Σ_i κ_i ∫ ds s^{-2} e^{-s m_i²}  ∫√g R
e diferem em como a integral em s (divergente em s → 0) é tratada.

1. Corte de tempo próprio único (Sakharov/Visser):
       1/G = Λ² Σκ_i / (2π)                         (já em induced_gravity.py)

2. Pauli-Villars (Sakharov finito, Adler-Zee): cada campo físico ganha reguladores de massa M_j e
   peso c_j com as regras de soma  Σ_j c_j = 0  e  Σ_j c_j M_j² = 0  (incluindo o campo físico).
   Então a integral é finita e
       1/G = (1/2π) Σ_i κ_i Σ_j c_ij M_ij² ln M_ij²
   Para dois reguladores (M1, M2) de um campo sem massa:
       1/G_i = (κ_i/2π) · M1² M2² ln(M1²/M2²) / (M1² - M2²)   → κ_i M²/(2π) quando M1, M2 → M.

3. Regularização dimensional / zeta: não há divergências de potência. Campos sem massa induzem
   1/G = 0; campos com massa dão um termo ∝ m² ln(m²/μ²), dependente da escala μ, que precisa de
   um contratermo 1/G_0: G deixa de ser predição.
       1/G(μ) = 1/G_0(μ) - (1/2π) Σ_i κ_i m_i² ln(m_i²/μ²)    (convenção de sinal da parte finita)
"""
from __future__ import annotations

from dataclasses import dataclass
from math import log, pi, sqrt

from .induced_gravity import KAPPA_VECTOR, KAPPA_WEYL, FieldContent, kappa_scalar


# ---------------------------------------------------------------- Pauli-Villars

def pv_weights(M1: float, M2: float) -> tuple[float, float]:
    """Pesos (c1, c2) de dois reguladores que, com o campo físico sem massa (c0 = 1), satisfazem
    1 + c1 + c2 = 0 e c1 M1² + c2 M2² = 0."""
    if M1 == M2:
        raise ValueError("use massas distintas (o limite M1 → M2 é tratado em pv_inverse_G_single)")
    d = M1**2 - M2**2
    return M2**2 / d, -(M1**2) / d


def pv_inverse_G_single(kappa: float, M1: float, M2: float) -> float:
    """1/G induzido por UM campo sem massa de peso κ, regulado por Pauli-Villars (M1, M2)."""
    if abs(M1 - M2) < 1e-12 * max(M1, M2):
        return kappa * M1**2 / (2 * pi)
    c1, c2 = pv_weights(M1, M2)
    return kappa * (c1 * M1**2 * log(M1**2) + c2 * M2**2 * log(M2**2)) / (2 * pi)


@dataclass(frozen=True)
class SpeciesRegulators:
    """Escala de regulador efetiva por tipo de campo (M_escalar, M_fermion, M_vetor)."""
    scalar: float = 1.0
    fermion: float = 1.0
    vector: float = 1.0


def pv_inverse_G(content: FieldContent, regs: SpeciesRegulators, split: float = 1.5) -> float:
    """1/G com Pauli-Villars, cada tipo de campo com seus reguladores (M, split·M)."""
    def one(k, M):
        return pv_inverse_G_single(float(k), M, split * M)
    return (content.real_scalars * one(kappa_scalar(content.xi), regs.scalar)
            + content.weyl_fermions * one(KAPPA_WEYL, regs.fermion)
            + content.gauge_vectors * one(KAPPA_VECTOR, regs.vector))


def critical_vector_regulator_ratio(content: FieldContent) -> float:
    """Razão M_vetor / M_matéria acima da qual 1/G fica negativo (mesmo split para todos)."""
    matter = content.real_scalars * float(kappa_scalar(content.xi)) + content.weyl_fermions * float(KAPPA_WEYL)
    vectors = -content.gauge_vectors * float(KAPPA_VECTOR)
    if vectors == 0:
        return float("inf")
    return sqrt(matter / vectors)


# ---------------------------------------------------------------- regularização dimensional / zeta

def dimreg_inverse_G_massless(content: FieldContent) -> float:
    """Campos sem massa: nenhuma escala, nenhuma divergência de potência → 1/G induzido = 0."""
    return 0.0


# massas em GeV (PDG, valores redondos) e pesos κ por grau de liberdade massivo
KAPPA_DIRAC = 2 * float(KAPPA_WEYL)          # +1/3
KAPPA_MASSIVE_VECTOR = float(KAPPA_VECTOR) + 1 / 6   # Proca: vetor - 1 fantasma de Stueckelberg = -1/2

SM_MASSIVE = [
    # (nome, massa GeV, κ, multiplicidade)
    ("top", 172.6, KAPPA_DIRAC, 3),
    ("bottom", 4.18, KAPPA_DIRAC, 3),
    ("charm", 1.27, KAPPA_DIRAC, 3),
    ("tau", 1.777, KAPPA_DIRAC, 1),
    ("W", 80.37, KAPPA_MASSIVE_VECTOR, 2),
    ("Z", 91.19, KAPPA_MASSIVE_VECTOR, 1),
    ("Higgs", 125.2, 1 / 6, 1),
]

PLANCK_MASS_GEV = 1.22e19


def dimreg_inverse_G_SM(mu: float) -> float:
    """Parte finita, dependente de μ (GeV), de 1/G induzido pelos campos massivos do MP (GeV²)."""
    return -sum(mult * k * m**2 * log(m**2 / mu**2) for _, m, k, mult in SM_MASSIVE) / (2 * pi)
