"""Modelo de Schwinger na rede: QED em 1+1D (U(1) DINÂMICO + férmions), diagonalização exata.

Férmions de Kogut-Susskind (staggered) em N sítios, condições abertas, Jordan-Wigner. A lei de Gauss
L_n - L_{n-1} = Q_n elimina o campo elétrico das ligações (Martinez et al., Nature 534, 516, 2016):

    H = w Σ_{n=0}^{N-2} (σ⁺_n σ⁻_{n+1} + h.c.) + (m/2) Σ_n (-1)^n σᶻ_n + J Σ_{n=0}^{N-2} L_n²
    L_n = ε₀ + Σ_{k≤n} Q_k,     Q_k = (σᶻ_k + (-1)^k) / 2     (carga no sítio k)

Ponto central: o termo J Σ L_n² acopla σᶻ_k σᶻ_l para TODOS os pares (k, l), com coeficiente
J·(N-1-max(k,l))/2 que NÃO decai com |k - l| (é o potencial linear, confinante, de 1+1D).
O Hamiltoniano de qubits é não local; a pergunta é se a geometria emergente do emaranhamento é local.

"Célula física" = par de sítios (2j, 2j+1) (elétron + pósitron), como a célula unitária da rede staggered.
"""
from __future__ import annotations

from itertools import combinations

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.linalg import eigsh


def sector_basis(N: int) -> np.ndarray:
    """Estados com carga total zero: N/2 spins para cima (σᶻ = +1 ⇔ bit = 1)."""
    states = [sum(1 << q for q in c) for c in combinations(range(N), N // 2)]
    return np.array(sorted(states), dtype=np.int64)


def zz_couplings(N: int, J: float) -> np.ndarray:
    """Matriz C_kl dos termos J·σᶻ_k σᶻ_l vindos de J Σ L_n² (para mostrar a não localidade)."""
    C = np.zeros((N, N))
    for k in range(N):
        for l in range(k + 1, N):
            C[k, l] = C[l, k] = J * (N - 1 - max(k, l)) / 2
    return C


def hamiltonian(N: int, w=1.0, m=0.0, J=1.0, eps0=0.0):
    basis = sector_basis(N)
    index = {int(s): i for i, s in enumerate(basis)}
    bits = ((basis[:, None] >> np.arange(N)) & 1).astype(float)
    sz = 2 * bits - 1
    stag = (-1.0) ** np.arange(N)
    Q = (sz + stag) / 2
    Ln = eps0 + np.cumsum(Q, axis=1)[:, :-1]
    diag = 0.5 * m * (sz * stag).sum(1) + J * (Ln**2).sum(1)
    rows, cols, vals = list(range(len(basis))), list(range(len(basis))), list(diag)
    for i, s in enumerate(basis):
        s = int(s)
        for n in range(N - 1):
            b1, b2 = (s >> n) & 1, (s >> (n + 1)) & 1
            if b1 != b2:
                t = s ^ ((1 << n) | (1 << (n + 1)))
                rows.append(i)
                cols.append(index[t])
                vals.append(w)
    H = coo_matrix((vals, (rows, cols)), shape=(len(basis),) * 2).tocsr()
    return H, basis


def ground_state(N: int, **kw):
    H, basis = hamiltonian(N, **kw)
    w, V = eigsh(H, k=1, which="SA", tol=1e-12)
    psi = np.zeros(1 << N)
    psi[basis] = V[:, 0]
    return float(w[0]), psi / np.linalg.norm(psi), basis


def free_fermion_energy(N: int, w=1.0, m=0.0) -> float:
    """J = 0: H = w Σ (c†c + h.c.) + m Σ (-1)^n n_n - (m/2)Σ(-1)^n; meio preenchimento (JW em 1D aberto é exato)."""
    h = np.zeros((N, N))
    for n in range(N - 1):
        h[n, n + 1] = h[n + 1, n] = w
    stag = (-1.0) ** np.arange(N)
    h += np.diag(m * stag)
    eps = np.sort(np.linalg.eigvalsh(h))
    return float(eps[: N // 2].sum() - 0.5 * m * stag.sum())


def expectation_charge(psi: np.ndarray, N: int) -> np.ndarray:
    idx = np.nonzero(np.abs(psi) > 0)[0]
    p = np.abs(psi[idx]) ** 2
    sz = 2 * ((idx[:, None] >> np.arange(N)) & 1) - 1
    return (p[:, None] * (sz + (-1.0) ** np.arange(N)) / 2).sum(0)
