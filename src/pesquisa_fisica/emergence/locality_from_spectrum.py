"""Localidade a partir do espectro (Cotler, Penington & Ranard, "Locality from the spectrum", CMP 2019).

Pergunta de H2: se o espaço emerge, de onde vem a própria divisão em "sítios" vizinhos? Cotler et al.
mostram que o ESPECTRO de um Hamiltoniano quase sempre determina se existe uma fatoração do espaço de
Hilbert em subsistemas na qual ele é local — e que, se existe, ela é essencialmente única.

Teste computável (contagem de parâmetros): um ansatz local de cadeia com p parâmetros tem p ≪ 2^n níveis.
  - Espectro de um Hamiltoniano local: o ajuste (só por autovalores!) chega a erro ~0.
  - Espectro genérico (GUE) ou de Hamiltoniano não local: o ajuste falha, com erro O(1).
Nada além da lista ordenada de autovalores entra no ajuste.

Ansatz: cadeia aberta de n qubits, invariante por translação, 2-local entre vizinhos:
    H(θ) = Σ_{i} Σ_{a,b∈{x,y,z}} J_ab σ^a_i σ^b_{i+1} + Σ_i Σ_a h_a σ^a_i      (9 + 3 = 12 parâmetros)
"""
from __future__ import annotations

from functools import reduce

import numpy as np
from scipy.optimize import least_squares

PAULI = {
    "x": np.array([[0, 1], [1, 0]], dtype=complex),
    "y": np.array([[0, -1j], [1j, 0]], dtype=complex),
    "z": np.array([[1, 0], [0, -1]], dtype=complex),
}
I2 = np.eye(2, dtype=complex)


def op_at(n: int, ops: dict[int, np.ndarray]) -> np.ndarray:
    return reduce(np.kron, [ops.get(i, I2) for i in range(n)])


def ansatz_terms(n: int) -> list[np.ndarray]:
    """Operadores (Hermitianos) multiplicados por cada parâmetro do ansatz TI 2-local."""
    terms = []
    for a in "xyz":
        for b in "xyz":
            terms.append(sum(op_at(n, {i: PAULI[a], i + 1: PAULI[b]}) for i in range(n - 1)))
    for a in "xyz":
        terms.append(sum(op_at(n, {i: PAULI[a]}) for i in range(n)))
    return terms


def hamiltonian(theta, terms) -> np.ndarray:
    return sum(t * T for t, T in zip(theta, terms))


def normalize(E: np.ndarray) -> np.ndarray:
    """Remove o deslocamento e a escala (só a FORMA do espectro importa)."""
    E = np.sort(np.real(E))
    E = E - E.mean()
    return E / np.linalg.norm(E)


def fit_spectrum(target: np.ndarray, terms, n_starts: int = 12, seed: int = 0) -> tuple[float, np.ndarray]:
    """Ajusta θ para que o espectro ordenado de H(θ) coincida com `target` (normalizado).
    Jacobiano por Hellmann-Feynman: ∂E_k/∂θ_j = ⟨ψ_k|T_j|ψ_k⟩. Devolve (erro RMS, θ)."""
    target = normalize(target)
    rng = np.random.default_rng(seed)

    def residual_and_jac(theta):
        H = hamiltonian(theta, terms)
        E, V = np.linalg.eigh(H)
        mu, nrm = E.mean(), np.linalg.norm(E - E.mean())
        En = (E - mu) / nrm
        # derivadas dos autovalores brutos
        dE = np.stack([np.real(np.einsum("ik,ij,jk->k", V.conj(), T, V)) for T in terms], axis=1)
        dE_c = dE - dE.mean(0)
        dn = (dE_c - np.outer(En, En @ dE_c)) / nrm               # derivada do espectro normalizado
        return En - target, dn

    best = (np.inf, None)
    for _ in range(n_starts):
        th0 = rng.standard_normal(len(terms))
        sol = least_squares(lambda t: residual_and_jac(t)[0], th0,
                            jac=lambda t: residual_and_jac(t)[1], xtol=1e-14, ftol=1e-14, gtol=1e-14, max_nfev=3000)
        err = float(np.sqrt(np.mean(sol.fun**2)))
        if err < best[0]:
            best = (err, sol.x)
        if best[0] < 1e-9:
            break
    return best


def gue_spectrum(dim: int, rng) -> np.ndarray:
    A = rng.standard_normal((dim, dim)) + 1j * rng.standard_normal((dim, dim))
    return np.linalg.eigvalsh((A + A.conj().T) / 2)


def random_nonlocal_hamiltonian_spectrum(n: int, rng, k_body: int = 3, n_terms: int = 12) -> np.ndarray:
    """Soma de n_terms strings de Pauli aleatórias de peso k_body em sítios ALEATÓRIOS (não vizinhos)."""
    H = np.zeros((2**n, 2**n), dtype=complex)
    for _ in range(n_terms):
        sites = rng.choice(n, size=k_body, replace=False)
        H += rng.standard_normal() * op_at(n, {int(s): PAULI[rng.choice(list("xyz"))] for s in sites})
    return np.linalg.eigvalsh(H)


def ansatz_terms_general(n: int) -> list[np.ndarray]:
    """Ansatz 2-local de vizinhos SEM invariância de translação: 9(n-1) + 3n parâmetros."""
    terms = []
    for i in range(n - 1):
        for a in "xyz":
            for b in "xyz":
                terms.append(op_at(n, {i: PAULI[a], i + 1: PAULI[b]}))
    for i in range(n):
        for a in "xyz":
            terms.append(op_at(n, {i: PAULI[a]}))
    return terms


def ansatz_terms_xxz_disordered(n: int) -> list[np.ndarray]:
    """Ansatz local SEM invariância de translação, enxuto: XX, YY, ZZ por ligação + campo z por sítio.
    3(n−1) + n parâmetros (n = 7: 25 parâmetros para 128 níveis)."""
    terms = []
    for i in range(n - 1):
        for a in "xyz":
            terms.append(op_at(n, {i: PAULI[a], i + 1: PAULI[a]}))
    for i in range(n):
        terms.append(op_at(n, {i: PAULI["z"]}))
    return terms


def random_nonlocal_sz_conserving_spectrum(n: int, rng, n_terms: int = 20) -> np.ndarray:
    """Controle não local que conserva S_z (mesma simetria do ansatz XXZ): trocas XX+YY e ZZ entre pares
    de sítios ALEATÓRIOS (não vizinhos) mais termos ZZZ de 3 corpos."""
    H = np.zeros((2**n, 2**n), dtype=complex)
    for _ in range(n_terms):
        i, j, k = rng.choice(n, size=3, replace=False)
        H += rng.standard_normal() * (op_at(n, {int(i): PAULI["x"], int(j): PAULI["x"]})
                                      + op_at(n, {int(i): PAULI["y"], int(j): PAULI["y"]}))
        H += rng.standard_normal() * op_at(n, {int(i): PAULI["z"], int(j): PAULI["z"], int(k): PAULI["z"]})
    return np.linalg.eigvalsh(H)
