"""Quench EXATO no Schwinger (setor de carga zero) e IM entre blocos arbitrários (partição de spin).

Base: bit n = sítio n, bit 1 = σᶻ = +1 (ocupado). Vácuo nu: sítios ímpares ocupados.
"""
from __future__ import annotations

import numpy as np
from scipy.sparse.linalg import eigsh, expm_multiply

from .schwinger import hamiltonian


def bare_vacuum_index(N: int) -> int:
    return sum(1 << k for k in range(1, N, 2))


def gap(N: int, m: float, J: float, w: float = 1.0) -> float:
    H, _ = hamiltonian(N, w=w, m=m, J=J)
    e = np.sort(eigsh(H, k=2, which="SA", tol=1e-10)[0])
    return float(e[1] - e[0])


def evolve(N: int, m: float, J: float, times, w: float = 1.0):
    """Lista de vetores de estado completos (2^N) nos tempos dados, a partir do vácuo nu."""
    H, basis = hamiltonian(N, w=w, m=m, J=J)
    v = np.zeros(len(basis), complex)
    v[np.searchsorted(basis, bare_vacuum_index(N))] = 1.0
    out, t_now = [], 0.0
    for t in times:
        if t > t_now:
            v = expm_multiply(-1j * (t - t_now) * H, v)
            t_now = t
        full = np.zeros(1 << N, complex)
        full[basis] = v
        out.append(full)
    return out


def _S(rho):
    w = np.linalg.eigvalsh(rho)
    w = w[w > 1e-14]
    return float(-(w * np.log(w)).sum())


def rdm(full: np.ndarray, N: int, sites) -> np.ndarray:
    """ρ dos sítios dados. Eixo k do tensor (2,)*N corresponde ao sítio N−1−k."""
    T = full.reshape((2,) * N)
    axes = [N - 1 - s for s in sites]
    rest = [k for k in range(N) if k not in axes]
    X = T.transpose(axes + rest).reshape(1 << len(sites), -1)
    return X @ X.conj().T


def block_mi(full, N, A, B) -> float:
    return _S(rdm(full, N, A)) + _S(rdm(full, N, B)) - _S(rdm(full, N, list(A) + list(B)))


def symmetric_blocks(N: int, lb: int, r: int):
    """A e B de tamanho lb, simétricos em torno do centro, separados por r sítios."""
    c = N // 2
    A = list(range(c - r // 2 - lb, c - r // 2))
    B = list(range(c + (r - r // 2), c + (r - r // 2) + lb))
    return A, B


def peak(times, I, thr=1e-3, drop=0.7):
    """(tem_pico, t_pico interpolado, I_max). Pico = máximo interior, I_max > thr, I(final) ≤ drop·I_max."""
    I = np.asarray(I); times = np.asarray(times)
    k = int(np.argmax(I))
    ok = 0 < k < len(I) - 1 and I[k] > thr and I[-1] <= drop * I[k]
    tp = times[k]
    if 0 < k < len(I) - 1:
        y0, y1, y2 = I[k - 1], I[k], I[k + 1]
        den = y0 - 2 * y1 + y2
        if den != 0:
            tp = times[k] + 0.5 * (y0 - y2) / den * (times[1] - times[0])
    return bool(ok), float(tp), float(I[k])


def classical_block_mi(full, N, A, B) -> float:
    """IM da distribuição de CARGAS (diagonal de ρ na base σᶻ): só observáveis locais invariantes de gauge."""
    H = lambda r: (lambda p: float(-(p[p > 1e-15] * np.log(p[p > 1e-15])).sum()))(np.real(np.diag(r)))
    return H(rdm(full, N, A)) + H(rdm(full, N, B)) - H(rdm(full, N, list(A) + list(B)))


def arrival_time(times, I, thr):
    I = np.asarray(I); times = np.asarray(times)
    k = np.where(I > thr)[0]
    if len(k) == 0 or k[0] == 0:
        return float("nan")
    k = k[0]
    return float(times[k - 1] + (thr - I[k - 1]) / (I[k] - I[k - 1]) * (times[k] - times[k - 1]))


# ---------------------------------------------------------------- controle: interação LOCAL sem campo de gauge

def local_interaction_hamiltonian(N: int, w: float = 1.0, m: float = 0.0, V: float = 0.0):
    """H_V = w Σ(σ⁺σ⁻ + h.c.) + (m/2) Σ (−1)^n σᶻ_n + V Σ n_i n_{i+1}  (Thirring/XXZ massivo na rede)."""
    from scipy.sparse import coo_matrix
    from .schwinger import sector_basis

    basis = sector_basis(N)
    index = {int(s): i for i, s in enumerate(basis)}
    bits = ((basis[:, None] >> np.arange(N)) & 1).astype(float)
    sz = 2 * bits - 1
    stag = (-1.0) ** np.arange(N)
    diag = 0.5 * m * (sz * stag).sum(1) + V * (bits[:, :-1] * bits[:, 1:]).sum(1)
    rows, cols, vals = list(range(len(basis))), list(range(len(basis))), list(diag)
    for i, s in enumerate(basis):
        s = int(s)
        for n in range(N - 1):
            if ((s >> n) & 1) != ((s >> (n + 1)) & 1):
                rows.append(i); cols.append(index[s ^ ((1 << n) | (1 << (n + 1)))]); vals.append(w)
    return coo_matrix((vals, (rows, cols)), shape=(len(basis),) * 2).tocsr(), basis


def gap_of(H) -> float:
    e = np.sort(eigsh(H, k=2, which="SA", tol=1e-10)[0])
    return float(e[1] - e[0])


def evolve_H(H, basis, N: int, times):
    v = np.zeros(len(basis), complex)
    v[np.searchsorted(basis, bare_vacuum_index(N))] = 1.0
    out, t_now = [], 0.0
    for t in times:
        if t > t_now:
            v = expm_multiply(-1j * (t - t_now) * H, v)
            t_now = t
        full = np.zeros(1 << N, complex)
        full[basis] = v
        out.append(full)
    return out


def truncated_linear_hamiltonian(N: int, w: float = 1.0, m: float = 0.0, J: float = 1.0, R: int = 1):
    """H_R = hopping + (m/2)Σ(−1)ⁿσᶻ − (J/2) Σ_{k,l} Q_k Q_l min(|k−l|, R), sem campo de gauge.
    Identidade (Q_tot = 0): J Σ_n L_n² = −(J/2) Σ Q_k Q_l |k−l| ⇒ R ≥ N−1 reproduz o Schwinger."""
    from scipy.sparse import diags
    H0, basis = local_interaction_hamiltonian(N, w, m, 0.0)
    bits = ((basis[:, None] >> np.arange(N)) & 1).astype(float)
    Q = bits - (np.arange(N) % 2)
    k = np.arange(N)
    W = np.minimum(np.abs(k[:, None] - k[None, :]), R).astype(float)
    pot = -0.5 * J * np.einsum("si,ij,sj->s", Q, W, Q)
    return (H0 + diags(pot)).tocsr(), basis
