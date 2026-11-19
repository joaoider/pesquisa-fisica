"""Gravidade induzida de Sakharov com corte de tempo próprio.

Ação efetiva euclidiana de um laço (bósons):
    Γ = -(1/2) ∫_{1/Λ²}^∞ (ds/s) Tr e^{-sD}
Com Tr e^{-sD} ≈ (4πs)^{-2} ∫√g Σ_k s^k a_k, o termo a_1 dá
    Γ ⊃ -(Λ²/32π²) ∫√g tr a_1.
Comparando com a ação de Einstein-Hilbert euclidiana S_E = -(1/16πG) ∫√g (R - 2Λ_cc):
    1/(16πG) = (Λ²/32π²) Σ_i κ_i    ⇒    1/G = Λ² Σ_i κ_i / (2π)
em que κ_i é o coeficiente de R em tr a_1 do campo i, já com o sinal estatístico
(férmions entram com sinal oposto) e com os fantasmas de Faddeev-Popov dos vetores.

Coeficientes κ por campo (unidades naturais ħ = c = 1):
    escalar real, acoplamento ξ:  1/6 - ξ
    férmion de Weyl:              +1/6      (Dirac: +1/3; tr a_1 = 4(1/6 - 1/4) R, sinal -1, metade p/ Weyl)
    vetor de gauge sem massa:     -2/3      (vetor: 4/6 - 1 = -1/3; menos 2 fantasmas escalares: -2/6)

Em "unidades de escalar" (κ / (1/6)) isso é N_0 + N_{1/2} - 4 N_1, a estrutura citada por Visser (2002).

ATENÇÃO: o resultado depende do esquema de regularização (um único Λ para todas as espécies,
corte de tempo próprio). Outros esquemas mudam os pesos relativos. Isso é parte do que se testa.
"""
from __future__ import annotations

from dataclasses import dataclass
from fractions import Fraction
from math import pi, sqrt

KAPPA_WEYL = Fraction(1, 6)
KAPPA_VECTOR = Fraction(-2, 3)


def kappa_scalar(xi: Fraction | float = 0) -> Fraction | float:
    return Fraction(1, 6) - xi


@dataclass(frozen=True)
class FieldContent:
    real_scalars: int = 0
    weyl_fermions: int = 0
    gauge_vectors: int = 0
    xi: Fraction | float = 0  # acoplamento não mínimo dos escalares

    def kappa_total(self):
        return (self.real_scalars * kappa_scalar(self.xi)
                + self.weyl_fermions * KAPPA_WEYL
                + self.gauge_vectors * KAPPA_VECTOR)

    def critical_xi(self) -> Fraction | None:
        """ξ acima do qual 1/G fica negativo (gravidade repulsiva)."""
        if self.real_scalars == 0:
            return None
        rest = self.weyl_fermions * KAPPA_WEYL + self.gauge_vectors * KAPPA_VECTOR
        return Fraction(1, 6) + Fraction(rest) / self.real_scalars


def standard_model(xi: Fraction | float = 0, right_handed_neutrinos: bool = False) -> FieldContent:
    """Conteúdo de campos do Modelo Padrão.

    Escalares: dubleto de Higgs = 4 escalares reais.
    Férmions de Weyl por geração: Q (3 cores × 2) + u^c (3) + d^c (3) + L (2) + e^c (1) = 15 (+1 com ν_R).
    Vetores: SU(3) 8 + SU(2) 3 + U(1) 1 = 12.
    """
    per_gen = 16 if right_handed_neutrinos else 15
    return FieldContent(real_scalars=4, weyl_fermions=3 * per_gen, gauge_vectors=12, xi=xi)


def inverse_G(content: FieldContent, cutoff: float) -> float:
    """1/G induzido (unidades naturais) para corte Λ = cutoff."""
    return cutoff**2 * float(content.kappa_total()) / (2 * pi)


def cutoff_for_newton(content: FieldContent, planck_mass: float = 1.0) -> float:
    """Λ necessário para reproduzir G = 1/M_Pl². Erro se 1/G induzido ≤ 0."""
    k = float(content.kappa_total())
    if k <= 0:
        raise ValueError(f"Σκ = {k:.4f} ≤ 0: o conteúdo de campos induz 1/G ≤ 0")
    return planck_mass * sqrt(2 * pi / k)


def induced_vacuum_energy_ratio(cutoff_over_planck: float, n_dof_signed: float) -> float:
    """|ρ_vac induzida| / ρ_Λ observada, com ρ_vac ~ N Λ⁴ / (64π²) (termo a_0).

    ρ_Λ observada ≈ 2.5e-47 GeV⁴ ≈ 7e-123 M_Pl⁴ (M_Pl = 1.22e19 GeV).
    """
    rho_obs_planck = 2.5e-47 / (1.22e19) ** 4
    rho_ind = abs(n_dof_signed) * cutoff_over_planck**4 / (64 * pi**2)
    return rho_ind / rho_obs_planck
