"""Cone de luz e assinatura lorentziana emergentes (H2).

1. Lieb-Robinson (1972): em Hamiltonianos locais, ‖[A_0(t), B_r]‖ é exponencialmente pequeno fora de
   r ≲ v t. Isso define uma ESTRUTURA CAUSAL: tempo e espaço entram com sinais opostos em s² = (vt)² − r².
   - férmions livres: {c_0(t), c_r†} = U_{0r}(t), U = e^{-iht} (exato, cadeias longas);
   - cadeia de Heisenberg (interagente), diagonalização exata: ‖[Z_0(t), Z_r]‖_F.
2. Isotropia de Lorentz emergente: perto do ponto de Dirac do fluxo π, E(q) ≈ v|q| é isotrópico
   (a rede quebra Lorentz só em ordem q²); um metal (Φ = 0) com superfície de Fermi não tem isso.
"""
from __future__ import annotations

from functools import reduce

import numpy as np
from scipy.linalg import expm


# ---------------------------------------------------------------- férmions livres 1D

def chain_hopping(N: int, w: float = 1.0) -> np.ndarray:
    h = np.zeros((N, N))
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = -w
    return h


def propagator_front(N: int, times, w: float = 1.0, origin: int | None = None) -> np.ndarray:
    """|U_{0,r}(t)| para r = 0..N-1-origin (origem no centro por padrão). Forma (len(times), n_r)."""
    origin = N // 2 if origin is None else origin
    h = chain_hopping(N, w)
    e, V = np.linalg.eigh(h)
    out = []
    for t in times:
        U = (V * np.exp(-1j * e * t)) @ V.T
        out.append(np.abs(U[origin, origin:]))
    return np.array(out)


def arrival_times(amplitudes: np.ndarray, times, threshold: float) -> np.ndarray:
    """Primeiro instante em que |U_{0r}| passa do limiar, para cada r (nan se nunca)."""
    arr = np.full(amplitudes.shape[1], np.nan)
    for r in range(amplitudes.shape[1]):
        idx = np.nonzero(amplitudes[:, r] > threshold)[0]
        if len(idx):
            arr[r] = times[idx[0]]
    return arr


# ---------------------------------------------------------------- cadeia de Heisenberg (interagente)

_X = np.array([[0, 1], [1, 0]], dtype=complex)
_Y = np.array([[0, -1j], [1j, 0]], dtype=complex)
_Z = np.array([[1, 0], [0, -1]], dtype=complex)


def _site_op(n, i, op):
    return reduce(np.kron, [op if k == i else np.eye(2) for k in range(n)])


def heisenberg(n: int, J: float = 1.0, Delta: float = 1.0) -> np.ndarray:
    H = np.zeros((2**n, 2**n), dtype=complex)
    for i in range(n - 1):
        H += J * (_site_op(n, i, _X) @ _site_op(n, i + 1, _X) + _site_op(n, i, _Y) @ _site_op(n, i + 1, _Y)
                  + Delta * _site_op(n, i, _Z) @ _site_op(n, i + 1, _Z))
    return H


def commutator_norms(n: int, times, H=None) -> np.ndarray:
    """C(r, t) = ‖[Z_0(t), Z_r]‖_F / ‖Z_0 Z_r‖_F (norma de Frobenius normalizada), origem no sítio 0."""
    H = heisenberg(n) if H is None else H
    e, V = np.linalg.eigh(H)
    Z0 = _site_op(n, 0, _Z)
    Zr = [_site_op(n, r, _Z) for r in range(n)]
    norm = np.sqrt(2**n)
    out = []
    for t in times:
        U = (V * np.exp(-1j * e * t)) @ V.conj().T
        Z0t = U.conj().T @ Z0 @ U
        out.append([np.linalg.norm(Z0t @ Zr[r] - Zr[r] @ Z0t) / (2 * norm) for r in range(n)])
    return np.array(out)


# ---------------------------------------------------------------- isotropia de Lorentz emergente 2D

def pi_flux_dispersion(kx, ky, t: float = 1.0):
    """Banda superior do fluxo π na rede quadrada: E = 2t √(cos²kx + cos²ky); Dirac em (π/2, π/2)."""
    return 2 * t * np.sqrt(np.cos(kx) ** 2 + np.cos(ky) ** 2)


def zero_flux_dispersion(kx, ky, t: float = 1.0):
    return -2 * t * (np.cos(kx) + np.cos(ky))


def dirac_anisotropy(q: float, n_dir: int = 64, t: float = 1.0) -> float:
    """(max − min)/média de E em torno do ponto de Dirac, num círculo de raio q."""
    ang = np.linspace(0, 2 * np.pi, n_dir, endpoint=False)
    E = pi_flux_dispersion(np.pi / 2 + q * np.cos(ang), np.pi / 2 + q * np.sin(ang), t)
    return float((E.max() - E.min()) / E.mean())


def fermi_velocity_anisotropy_zero_flux(n_dir: int = 256, t: float = 1.0, mu: float = -1.0) -> float:
    """Metal (Φ = 0): anisotropia da velocidade de Fermi |∇E| ao longo da superfície de Fermi E = μ."""
    ang = np.linspace(0, 2 * np.pi, n_dir, endpoint=False)
    vs = []
    for a in ang:
        # raio k_F(a) por bissecção ao longo da direção a
        lo, hi = 0.0, np.pi
        for _ in range(60):
            mid = (lo + hi) / 2
            if zero_flux_dispersion(mid * np.cos(a), mid * np.sin(a), t) < mu:
                lo = mid
            else:
                hi = mid
        kx, ky = lo * np.cos(a), lo * np.sin(a)
        vs.append(2 * t * np.hypot(np.sin(kx), np.sin(ky)))
    vs = np.array(vs)
    return float((vs.max() - vs.min()) / vs.mean())
