"""Tempo emergente: mecanismo de Page-Wootters (1983).

O universo inteiro (relógio C + sistema S) está num estado ESTACIONÁRIO:
    (H_C + H_S) |Ψ> = 0       → nada "acontece" globalmente.
Mesmo assim, o estado do sistema condicionado à leitura t do relógio,
    |ψ(t)> ∝ <t|_C |Ψ>,
obedece à equação de Schrödinger: o tempo surge como correlação (emaranhamento) entre partes.

Modelo: relógio de N níveis com energias E_k = -k ω (k = 0..N-1) e estados de leitura
    |t> = (1/√N) Σ_k e^{-i E_k t} |k>;
sistema com espectro em múltiplos inteiros de ω.
"""
from __future__ import annotations

import numpy as np


def build(system_levels: np.ndarray, psi0: np.ndarray, n_clock: int, omega: float = 1.0):
    """Devolve (Ψ global no kernel de H_total, H_total, energias do relógio).

    system_levels: inteiros n_j (H_S diagonal com autovalores n_j ω), 0 ≤ n_j < n_clock.
    psi0: estado inicial do sistema na base de energia.
    """
    system_levels = np.asarray(system_levels, dtype=int)
    E_c = -omega * np.arange(n_clock)
    E_s = omega * system_levels
    H_tot = np.add.outer(E_c, E_s).ravel()                    # diagonal: C ⊗ S
    # Ψ = Σ_j ψ0_j |k = n_j> ⊗ |j>  (cada componente do sistema pareada com o nível do relógio que zera a energia)
    Psi = np.zeros((n_clock, len(system_levels)), dtype=complex)
    for j, n in enumerate(system_levels):
        Psi[n, j] = psi0[j]
    Psi = Psi / np.linalg.norm(Psi)
    return Psi, H_tot, E_c


def clock_state(t: float, E_c: np.ndarray) -> np.ndarray:
    return np.exp(-1j * E_c * t) / np.sqrt(len(E_c))


def conditional_state(Psi: np.ndarray, t: float, E_c: np.ndarray) -> np.ndarray:
    phi = clock_state(t, E_c).conj() @ Psi
    return phi / np.linalg.norm(phi)


def schrodinger(psi0: np.ndarray, system_levels, t: float, omega: float = 1.0) -> np.ndarray:
    return np.exp(-1j * omega * np.asarray(system_levels) * t) * psi0 / np.linalg.norm(psi0)


def fidelity(a: np.ndarray, b: np.ndarray) -> float:
    return float(abs(np.vdot(a, b)) ** 2)


def clock_system_entanglement(Psi: np.ndarray) -> float:
    """Entropia de emaranhamento relógio-sistema (Schmidt). Zero ⇒ não há tempo emergente."""
    s = np.linalg.svd(Psi, compute_uv=False) ** 2
    s = s[s > 1e-15]
    return float(-np.sum(s * np.log(s)))
