"""Rede tensorial holográfica com tensores perfeitos (Pastawski, Yoshida, Harlow & Preskill — "HaPPY", 2015).

Bloco básico: o código de 5 qubits [[5,1,3]] visto como um tensor de 6 pernas T[l, p1..p5]
(l = qubit lógico, p = qubits físicos). Ele é PERFEITO: qualquer divisão 3|3 das pernas dá uma matriz
proporcional a uma unitária ⇒ o tensor é isometria de qualquer conjunto menor para o complementar.

Rede (árvore de 2 camadas, "disco hiperbólico" mínimo):
    T0 (centro): pernas (b0, a1..a5); b0 = qubit de bulk central (mantido como referência)
    Ti, i = 1..5: pernas (a_i, b_i, s_i1..s_i4); b_i = bulk fixado em |0⟩; s = 20 qubits de BORDA
Estado: |Ψ⟩ em C²(b0) ⊗ C^{2^20}(borda), normalizado.

Testes de holografia:
  - Ryu-Takayanagi discreto: S(A) = (corte mínimo na rede) · ln 2 para regiões A da borda.
  - Reconstrução da cunha de emaranhamento: o qubit de bulk b0 é recuperável de A ⇔ I(b0:A) = 2 ln 2.
"""
from __future__ import annotations

from itertools import product

import numpy as np

X = np.array([[0, 1], [1, 0]], dtype=complex)
Z = np.array([[1, 0], [0, -1]], dtype=complex)
I2 = np.eye(2, dtype=complex)


def _kron(ops):
    out = np.array([[1.0 + 0j]])
    for o in ops:
        out = np.kron(out, o)
    return out


def five_qubit_code_tensor() -> np.ndarray:
    """T[l, p1, ..., p5] com l ∈ {0,1}: |l_L⟩ do código de 5 qubits (estabilizadores XZZXI cíclicos)."""
    gens = []
    base = ["X", "Z", "Z", "X", "I"]
    for s in range(4):
        word = base[-s:] + base[:-s] if s else base
        gens.append(_kron([{"X": X, "Z": Z, "I": I2}[c] for c in word]))
    P = np.eye(32, dtype=complex)
    for g in gens:
        P = P @ (np.eye(32) + g) / 2
    # |0_L⟩ ∝ P|00000⟩ ; |1_L⟩ = X⊗5 |0_L⟩ (X lógico)
    zero = P[:, 0]
    zero = zero / np.linalg.norm(zero)
    one = _kron([X] * 5) @ zero
    T = np.stack([zero, one]).reshape([2] * 6)
    return T / np.linalg.norm(T) * np.sqrt(2)          # normalização: Σ|T|² = 2 (isometria lógico→físico)


def is_perfect(T: np.ndarray, tol: float = 1e-10) -> bool:
    """Toda divisão 3|3 das 6 pernas dá matriz ∝ unitária."""
    from itertools import combinations
    for A in combinations(range(6), 3):
        B = [k for k in range(6) if k not in A]
        M = np.transpose(T, list(A) + B).reshape(8, 8)
        G = M.conj().T @ M
        if not np.allclose(G, G[0, 0] * np.eye(8), atol=tol):
            return False
    return True


def build_network_state() -> np.ndarray:
    """|Ψ⟩ com forma (2, 2^20): eixo 0 = b0, eixo 1 = 20 qubits de borda na ordem cíclica (i, j)."""
    T = five_qubit_code_tensor()
    # T0[b0, a1..a5]
    T0 = T
    # Ti com b_i = |0⟩: Ti[a, s1..s4] = T[b=0? ] — usamos a perna 0 (lógica) como b_i e a perna 1 como a_i
    Ti = T[0]                                   # fixa b_i = 0 → pernas (a_i, s1..s4)
    psi = np.einsum("zabcde,aFGHI,bJKLM,cNOPQ,dRSTU,eVWXY->zFGHIJKLMNOPQRSTUVWXY",
                    T0, Ti, Ti, Ti, Ti, Ti, optimize=True)
    psi = psi.reshape(2, 2**20)
    return psi / np.linalg.norm(psi)


# ---------------------------------------------------------------- entropias

def entropy_of_boundary_region(psi: np.ndarray, region: list[int], include_bulk: bool = False) -> float:
    """S(A) em unidades de ln 2, A ⊂ {0..19} (qubits de borda), opcionalmente junto com b0."""
    t = psi.reshape([2] * 21)                    # eixo 0 = b0, eixos 1..20 = borda
    keep = ([0] if include_bulk else []) + [1 + q for q in region]
    rest = [k for k in range(21) if k not in keep]
    M = np.transpose(t, keep + rest).reshape(2 ** len(keep), -1)
    if M.shape[0] > M.shape[1]:
        M = M.T
    p = np.linalg.svd(M, compute_uv=False) ** 2
    p = p[p > 1e-14]
    return float(-np.sum(p * np.log2(p)))


def min_cut(region: list[int]) -> int:
    """Corte mínimo (nº de pernas cortadas) separando A do resto (borda complementar + b0), por força bruta
    sobre a atribuição dos 6 tensores aos lados. Pernas b_i fixadas em |0⟩ não contam (não carregam estado)."""
    A = set(region)
    best = 10**9
    for sides in product([0, 1], repeat=6):           # 1 = lado de A; nós: T0, T1..T5
        cut = 0
        cut += 1 if sides[0] == 1 else 0              # perna b0 (pertence ao complemento)
        for i in range(1, 6):
            cut += sides[0] != sides[i]               # ligação a_i
            for j in range(4):
                q = 4 * (i - 1) + j
                in_A = q in A
                cut += (sides[i] == 1) != in_A        # perna de borda do lado errado
        best = min(best, cut)
    return best


def contiguous_regions(length: int, n: int = 20) -> list[list[int]]:
    return [[(s + k) % n for k in range(length)] for s in range(n)]
