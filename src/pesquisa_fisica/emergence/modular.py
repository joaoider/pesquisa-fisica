"""Hamiltoniano de emaranhamento (modular) e primeira lei do emaranhamento (H2 → Einstein).

Caminho de Jacobson (1995, 2015) e Faulkner-Guica-Hartman-Myers-Van Raamsdonk (2014): se
  (i)  o Hamiltoniano de emaranhamento K do vácuo restrito a uma região é LOCAL e GEOMÉTRICO
       (Bisognano-Wichmann: K = ∫ β(x) T₀₀(x) dx, com β(x) = 2π (x)(ℓ−x)/(ℓ v) num intervalo de
       comprimento ℓ: a "temperatura de Unruh" local do vácuo), e
  (ii) vale a primeira lei do emaranhamento δS = δ⟨K⟩ para perturbações pequenas,
então, num contexto holográfico, impor (ii) para todas as regiões equivale às equações de Einstein
linearizadas. Aqui testamos (i) e (ii) numa matéria de rede: cadeia crítica de férmions livres
(CFT com c = 1, v = 2w sin k_F = 2 a meio preenchimento).

Para estados gaussianos: ρ_A ∝ exp(−Σ h_ij c_i† c_j) com h = ln[(1 − C_A)/C_A] (C_A = ⟨c_i† c_j⟩ em A).
O logaritmo é mal condicionado (autovalores exponencialmente próximos de 0 e 1): usamos mpmath.
"""
from __future__ import annotations

import mpmath as mp
import numpy as np


def infinite_chain_correlation(ell: int, kF: float = np.pi / 2, dps: int = 60):
    """C_ij = sin(kF (i−j)) / (π (i−j)), C_ii = kF/π (cadeia infinita, exata), em precisão mp."""
    mp.mp.dps = dps
    kF = mp.mpf(kF)
    C = mp.matrix(ell, ell)
    for i in range(ell):
        for j in range(ell):
            C[i, j] = kF / mp.pi if i == j else mp.sin(kF * (i - j)) / (mp.pi * (i - j))
    return C


def entanglement_hamiltonian(C_mp) -> np.ndarray:
    """h = ln[(1 − C)/C] calculado em precisão estendida; devolve float64."""
    E, Q = mp.eighe(C_mp) if hasattr(mp, "eighe") else mp.eigsy(C_mp)
    n = C_mp.rows
    D = mp.matrix(n, n)
    for k in range(n):
        lam = E[k]
        D[k, k] = mp.log((1 - lam) / lam)
    H = Q * D * Q.T
    return np.array(H.tolist(), dtype=float)


def bw_profile(ell: int, v: float = 2.0) -> np.ndarray:
    """β nas ligações (x + ½), x = 0..ℓ−2, com o intervalo ocupando [−½, ℓ − ½]."""
    x = np.arange(ell - 1) + 1.0      # posição da ligação medida a partir da borda esquerda (−½)
    return 2 * np.pi * x * (ell - x) / (ell * v)


# ---------------------------------------------------------------- primeira lei em cadeia finita

def chain_correlation(N: int, potential=None, w: float = 1.0) -> np.ndarray:
    h = np.zeros((N, N))
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = -w
    if potential is not None:
        h += np.diag(potential)
    e, U = np.linalg.eigh(h)
    occ = U[:, : N // 2]
    return occ @ occ.T


def fermion_entropy(C: np.ndarray) -> float:
    lam = np.clip(np.linalg.eigvalsh(C), 1e-15, 1 - 1e-15)
    return float(-np.sum(lam * np.log(lam) + (1 - lam) * np.log(1 - lam)))


def first_law_check(N: int, ell: int, potential: np.ndarray, eps: float, v: float = 2.0):
    """Devolve (δS, δ⟨K_BW⟩) para a perturbação ε·V do Hamiltoniano, intervalo central de comprimento ℓ.
    K_BW = Σ_ligações β(x+½)·(−(c_x† c_{x+1} + h.c.))  (densidade de energia de ligação, sem termo de μ
    porque o meio preenchimento é simétrico partícula-buraco)."""
    a = (N - ell) // 2
    A = slice(a, a + ell)
    C0 = chain_correlation(N)
    C1 = chain_correlation(N, eps * potential)
    dS = fermion_entropy(C1[A, A]) - fermion_entropy(C0[A, A])
    beta = bw_profile(ell, v)
    bond = lambda C: np.array([-2 * C[a + x, a + x + 1] for x in range(ell - 1)])
    dK = float(np.sum(beta * (bond(C1) - bond(C0))))
    return dS, dK


def chain_correlation_bonds(N: int, bond_mod=None, w: float = 1.0) -> np.ndarray:
    """Cadeia com salto modulado w_x = w (1 + f_x): perturbação da 'métrica' local (velocidade local)."""
    h = np.zeros((N, N))
    f = np.zeros(N - 1) if bond_mod is None else bond_mod
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = -w * (1 + f[i])
    e, U = np.linalg.eigh(h)
    occ = U[:, : N // 2]
    return occ @ occ.T


def first_law_bonds(N: int, ell: int, f: np.ndarray, eps: float, K_exact: np.ndarray | None = None, v: float = 2.0):
    """(δS, δ⟨K_BW⟩, δ⟨K_exato⟩) para a modulação de salto ε·f, intervalo central de comprimento ℓ.
    K_exato: Hamiltoniano de emaranhamento de partícula única do vácuo (ℓ×ℓ), δ⟨K⟩ = Σ h_ij δC_ij."""
    a = (N - ell) // 2
    A = slice(a, a + ell)
    C0 = chain_correlation_bonds(N)
    C1 = chain_correlation_bonds(N, eps * f)
    dS = fermion_entropy(C1[A, A]) - fermion_entropy(C0[A, A])
    beta = bw_profile(ell, v)
    bond = lambda C: np.array([-2 * C[a + x, a + x + 1] for x in range(ell - 1)])
    dK_bw = float(np.sum(beta * (bond(C1) - bond(C0))))
    dK_ex = None if K_exact is None else float(np.sum(K_exact * (C1[A, A] - C0[A, A])))
    return dS, dK_bw, dK_ex


# ---------------------------------------------------------------- δS para todos os intervalos: "espaço cinemático" dS₂

def chain_correlation_general(N: int, hopping: np.ndarray) -> np.ndarray:
    h = np.zeros((N, N))
    for i in range(N - 1):
        h[i, i + 1] = h[i + 1, i] = -hopping[i]
    e, U = np.linalg.eigh(h)
    occ = U[:, : N // 2]
    return occ @ occ.T


def delta_S_grid(C0: np.ndarray, C1: np.ndarray, centers: np.ndarray, radii: np.ndarray) -> np.ndarray:
    """δS(x, R) para intervalos [x − R + ½, x + R − ½] (2R sítios, centro x numa ligação).
    centers: posições de ligação x (o intervalo começa no sítio x − R + 1 se x é a ligação entre x e x+1)."""
    out = np.zeros((len(centers), len(radii)))
    for i, x in enumerate(centers):
        for j, R in enumerate(radii):
            a = int(x) - int(R) + 1
            A = slice(a, a + 2 * int(R))
            out[i, j] = fermion_entropy(C1[A, A]) - fermion_entropy(C0[A, A])
    return out


def ds2_klein_gordon_residual(dS: np.ndarray, radii: np.ndarray) -> np.ndarray:
    """R²(−∂²_R + ∂²_x) δS + 2 δS por diferenças finitas (passo 1 em x e em R), pontos internos.
    Para uma CFT₂ com Bisognano-Wichmann (⇔ holograficamente: Einstein linearizado em AdS₃), é zero."""
    d2R = dS[1:-1, 2:] - 2 * dS[1:-1, 1:-1] + dS[1:-1, :-2]
    d2x = dS[2:, 1:-1] - 2 * dS[1:-1, 1:-1] + dS[:-2, 1:-1]
    R = radii[1:-1][None, :]
    return R**2 * (-d2R + d2x) + 2 * dS[1:-1, 1:-1]


def bw_kernel_matrix(centers, radii, y) -> np.ndarray:
    """δS(x,R) ≈ Σ_y K[(x,R), y] ε(y): núcleo β = 2π(R² − u²)/(2R v), v = 2 ⇒ π(R² − u²)/(2R),
    u = y − x (ε = variação da energia da ligação y)."""
    rows = []
    for x in centers:
        for R in radii:
            u = y - x
            rows.append(np.where(np.abs(u) < R, np.pi * (R**2 - u**2) / (R * 2.0), 0.0))
    return np.array(rows)


# ---------------------------------------------------------------- bulk AdS₃ reconstruído e geodésicas (RT)

def central_charge_from_chain(N: int = 2000, lengths=np.arange(20, 200, 10)) -> float:
    """c pelo ajuste S = (c/6) ln[(2N/π) sin(πℓ/N)] + const no meio da cadeia aberta (Calabrese-Cardy, 1 borda
    distante: usamos intervalos centrais e a forma de cadeia infinita S = (c/3) ln ℓ + const)."""
    C = chain_correlation_general(N, np.ones(N - 1))
    S = [fermion_entropy(C[N // 2 - L // 2: N // 2 - L // 2 + L, N // 2 - L // 2: N // 2 - L // 2 + L]) for L in lengths]
    return float(3 * np.polyfit(np.log(lengths), S, 1)[0])


def reconstruct_energy_density(dS_small: np.ndarray, centers, radii_small, y, reg: float = 1e-6) -> np.ndarray:
    """Inverte δS(x, R) = Σ_y K ε(y) usando SÓ intervalos pequenos (mínimos quadrados com Tikhonov)."""
    K = bw_kernel_matrix(centers, radii_small, y)            # já inclui 1/v (v = 2)
    b = dS_small.ravel()
    A = K.T @ K + reg * np.trace(K.T @ K) / K.shape[1] * np.eye(K.shape[1])
    return np.linalg.solve(A, K.T @ b)


def rt_geodesic_delta_S(eps_y: np.ndarray, y: np.ndarray, x: float, R: float, G: float, n_theta: int = 4000) -> float:
    """δS = δL/(4G) pela variação do comprimento da geodésica semicircular z² + (y−x)² = R² no AdS₃ com
    métrica perturbada g = [dz² + (1 + z² t_xx) dx² + …]/z², t_xx(y) = 8πG ε(y)/v  (Fefferman-Graham, ℓ_AdS = 1).
    Integral feita numericamente em θ, independente da fórmula do núcleo."""
    th = np.linspace(1e-6, np.pi - 1e-6, n_theta)
    yy = x + R * np.cos(th)
    t_xx = 8 * np.pi * G * np.interp(yy, y, eps_y, left=0.0, right=0.0) / 2.0
    # δL = ½ ∫ δg_xx (dx/dθ)² / |dX/dθ|_g dθ, δg_xx = t_xx (o fator z² cancela com 1/z²)
    integrand = 0.5 * t_xx * (R * np.sin(th)) ** 2 * np.sin(th)
    dL = np.trapezoid(integrand, th)
    return float(dL / (4 * G))
