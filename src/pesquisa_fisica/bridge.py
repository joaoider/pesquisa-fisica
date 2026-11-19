"""Ponte H1 ↔ H2: a entropia de emaranhamento de área é S = A/(4G) com o MESMO G induzido de H1?

Tudo no mesmo esquema: corte de tempo próprio (s > 1/Λ²), como em induced_gravity.py.

Entropia de emaranhamento através de um plano de área A (método das réplicas): ln Z_α no espaço C_α × R²,
em que C_α é um cone de ângulo 2πα. O heat kernel de um escalar no cone tem um termo extra da ponta,
    Tr K_{C_α} ⊃ (1/12)(1/α − α),
e S = (1 − α∂_α) ln Z_α |_{α=1} dá, por escalar mínimo,
    S_escalar = A Λ² / (48π).
Em H1: 1/G = Λ² Σκ/(2π) com κ_escalar = 1/6 ⇒ A/(4G) = A Λ²/(48π) por escalar. Iguais.

Pesos por espécie em "unidades de escalar":
                       entropia (contagem de estados)   A/(4G_ind)          diferença
  escalar ξ = 0                  1                       1                   0
  escalar ξ ≠ 0                  1                       1 − 6ξ             termo de Wald (−6ξ)
  férmion de Weyl                1                       1                   0
  vetor de gauge                 2 (polarizações)        −4                  termo de contato de Kabat (−6)
O termo de contato dos vetores é a contribuição dos modos de borda (edge modes) do campo de gauge
(Kabat 1995; Donnelly & Wall 2015/2016). Desde a v0.12 ele é CALCULADO aqui (seção 4: espectro exato da
1-forma no cone, com a holonomia da ponta), não tomado da literatura: −6,000 por vetor.
"""
from __future__ import annotations

from math import pi

import numpy as np
from scipy.special import jn_zeros

# ---------------------------------------------------------------- 1. heat kernel no cone, pelo espectro exato


def cone_heat_trace(alpha_inv: int, s: np.ndarray, j_max: float | None = None) -> np.ndarray:
    """Tr e^{sΔ} (Dirichlet) num disco cônico de raio 1 e ângulo 2π/n (α = 1/n, n inteiro):
    autofunções J_{n|m|}(j r) e^{imφ}, autovalores j² com j zero de J_{n|m|}."""
    s = np.atleast_1d(np.asarray(s, dtype=float))
    j_max = j_max or np.sqrt(40 / s.min())
    total = np.zeros_like(s)
    m = 0
    while True:
        order = alpha_inv * m
        if order > j_max:
            break
        nt = int(j_max / np.pi) + 5
        z = jn_zeros(order, nt)
        z = z[z <= j_max + 10]
        deg = 1 if m == 0 else 2
        total += deg * np.exp(-np.outer(s, z**2)).sum(1)
        m += 1
    return total


def fit_cone_constant(alpha_inv: int, s_min: float = 0.002, s_max: float = 0.02, n_pts: int = 40) -> dict:
    """Ajusta T(s) = a/s + b/√s + c + d√s + e s e devolve os coeficientes."""
    s = np.linspace(s_min, s_max, n_pts)
    T = cone_heat_trace(alpha_inv, s)
    V = np.stack([1 / s, 1 / np.sqrt(s), np.ones_like(s), np.sqrt(s), s], axis=1)
    a, b, c, d, e = np.linalg.lstsq(V, T, rcond=None)[0]
    return {"area": a, "perimetro": b, "constante": c}


def expected_cone_constant(alpha: float) -> float:
    """Constante esperada: borda circular α/6 (curvatura geodésica) + ponta (1/12)(1/α − α)."""
    return alpha / 6 + (1 / 12) * (1 / alpha - alpha)


# ---------------------------------------------------------------- 2. réplicas: entropia por escalar

def entanglement_entropy_per_area_scalar(cutoff: float) -> float:
    """S/A = Λ²/(48π): ln Z_α ⊃ ½∫_{1/Λ²} ds/s · (1/12)(1/α−α) · A/(4πs) ⇒ (1 − α∂_α) dá fator 2."""
    conical = lambda a: (1 / 12) * (1 / a - a)
    d_alpha = 1e-6
    dconical = (conical(1 + d_alpha) - conical(1 - d_alpha)) / (2 * d_alpha)
    replica_factor = conical(1.0) - 1.0 * dconical                   # (1 − α∂_α) em α = 1  → 1/6
    integral = 0.5 * (1 / (4 * pi)) * cutoff**2                      # ½ ∫_{1/Λ²}^∞ ds/s · 1/s
    return replica_factor * integral


# ---------------------------------------------------------------- 3. contabilidade por espécie

ENTROPY_WEIGHT = {"escalar": 1.0, "weyl": 1.0, "vetor": 2.0}          # contagem de estados (bulk)
KABAT_CONTACT = {"escalar": 0.0, "weyl": 0.0, "vetor": -6.0}           # modos de borda; conferido em gauge_field_entropy_weights()


def g_weight(species: str, xi: float = 0.0) -> float:
    """Peso de A/(4G_ind) em unidades de escalar: κ_i/(1/6) (de induced_gravity.py)."""
    return {"escalar": 1.0 - 6.0 * xi, "weyl": 1.0, "vetor": -4.0}[species]


def standard_model_tally(xi: float = 0.0, edge_modes: bool = False, right_handed_neutrinos: bool = False) -> dict:
    n = {"escalar": 4, "weyl": 48 if right_handed_neutrinos else 45, "vetor": 12}
    S = sum(n[k] * (ENTROPY_WEIGHT[k] + (KABAT_CONTACT[k] if edge_modes else 0.0)) for k in n)
    G = sum(n[k] * g_weight(k, xi) for k in n)
    if edge_modes:
        # termo de Wald do ξ entra junto com os termos de contato (não é contagem de estados)
        S += n["escalar"] * (-6.0 * xi)
    return {"S_emaranhamento": S, "A_sobre_4G": G, "razao": S / G if G else float("inf")}


# ---------------------------------------------------------------- 4. termo de contato de Kabat pelo espectro exato

def _bessel_trace(orders, s, j_max):
    """Σ_ν Σ_k e^{−s j²_{ν,k}} para uma lista de ordens inteiras ν ≥ 0 (com multiplicidade)."""
    total = np.zeros_like(s)
    for nu in orders:
        if nu > j_max:
            continue
        nt = int(j_max / np.pi) + 5
        z = jn_zeros(int(nu), nt)
        z = z[z <= j_max + 10]
        total += np.exp(-np.outer(s, z**2)).sum(1)
    return total


def cone_one_form_heat_trace(alpha_inv: int, s: np.ndarray, j_max: float | None = None) -> np.ndarray:
    """Tr e^{sΔ₁} para 1-formas no disco cônico de ângulo 2π/n (Dirichlet nas componentes V_x ± iV_y).

    No cone, V_x + iV_y = e^{iφ}(V_r + iV_φ): com (V_r + iV_φ) ∝ e^{i m φ/α} (univalente no cone), as componentes
    complexas carregam índice angular ν = m/α ± 1 — a holonomia do transporte do vetor em torno da ponta.
    Para α = 1/n: ν = n m + 1 (inteiro), m ∈ ℤ (o sinal − é o conjugado). Fator 2: campo vetorial real.
    """
    s = np.atleast_1d(np.asarray(s, dtype=float))
    j_max = j_max or np.sqrt(40 / s.min())
    M = int(j_max / alpha_inv) + 3
    orders = []
    for m in range(-M, M + 1):
        orders.append(abs(alpha_inv * m + 1))      # componente V_x + iV_y
    # cada autofunção complexa da componente V_x + iV_y dá 2 modos REAIS do campo vetorial (f e i·f)
    return 2 * _bessel_trace(orders, s, j_max)


def fit_constant(trace_fn, alpha_inv: int, s_min=0.002, s_max=0.02, n_pts=40) -> float:
    s = np.linspace(s_min, s_max, n_pts)
    T = trace_fn(alpha_inv, s)
    V = np.stack([1 / s, 1 / np.sqrt(s), np.ones_like(s), np.sqrt(s), s], axis=1)
    return float(np.linalg.lstsq(V, T, rcond=None)[0][2])


def decompose_tip(constants: dict[int, float]) -> dict:
    """Ajusta C(α) = b·α + p·(1/α − α) + q·(1 − α) a valores em α = 1/n (b: borda; p, q: ponta)."""
    ns = sorted(constants)
    a = np.array([1 / n for n in ns])
    A = np.stack([a, 1 / a - a, 1 - a], axis=1)
    sol, res, *_ = np.linalg.lstsq(A, np.array([constants[n] for n in ns]), rcond=None)
    pred = A @ sol
    return {"b": sol[0], "p": sol[1], "q": sol[2], "max_resid": float(np.abs(pred - np.array([constants[n] for n in ns])).max())}


def entropy_weight_from_tip(p: float, q: float) -> float:
    """Peso de entropia em unidades de escalar: (1 − α∂_α)[p(1/α − α) + q(1 − α)] em α = 1, dividido pelo do
    escalar ((1 − α∂_α)(1/12)(1/α − α) = 1/6). (1 − α∂_α)(1/α − α) = 2 e (1 − α∂_α)(1 − α) = 1."""
    return (2 * p + q) / (1 / 6)


def gauge_field_entropy_weights() -> dict:
    """Pesos (unidades de escalar) obtidos do espectro exato no cone, para o campo de gauge em 4D:
    2 componentes transversas (escalares) + 1-forma no plano do cone − 2 fantasmas (escalares)."""
    Cs = {n: fit_constant(cone_heat_trace, n) for n in (1, 2, 3, 4, 5)}
    Cv = {n: fit_constant(cone_one_form_heat_trace, n) for n in (1, 2, 3, 4, 5)}
    ds, dv = decompose_tip(Cs), decompose_tip(Cv)
    w_s = entropy_weight_from_tip(ds["p"], ds["q"])
    w_v = entropy_weight_from_tip(dv["p"], dv["q"])
    total = 2 * w_s + w_v - 2 * w_s
    fisico = 2 * w_s + 2 * w_s - 2 * w_s           # se a 1-forma valesse 2 escalares (sem holonomia)
    return {"escalar": w_s, "um_forma_no_cone": w_v, "campo_de_gauge": total,
            "contagem_fisica": fisico, "termo_de_contato": total - fisico,
            "q_holonomia": dv["q"], "ajuste_escalar": ds, "ajuste_um_forma": dv}


# ---------------------------------------------------------------- 5. parte universal na rede: cadeia radial l = 0

def radial_chain_l(N: int, l: int) -> tuple[np.ndarray, np.ndarray]:
    """Srednicki (1993): onda parcial l de um escalar livre em 3+1D numa rede radial (a = 1, j = 1..N):
    H = ½ Σ_j [π_j² + (j+½)² (φ_j/j − φ_{j+1}/(j+1))² + l(l+1)/j² φ_j²]. Devolve (X, P) do vácuo."""
    K = np.zeros((N, N))
    for j in range(1, N):
        c = (j + 0.5) ** 2
        a, b = 1 / j, -1 / (j + 1)
        K[j - 1, j - 1] += c * a * a
        K[j, j] += c * b * b
        K[j - 1, j] += c * a * b
        K[j, j - 1] += c * a * b
    K[N - 1, N - 1] += (N + 0.5) ** 2 / N**2            # parede em N+1 (φ_{N+1} = 0)
    for j in range(1, N + 1):
        K[j - 1, j - 1] += l * (l + 1) / j**2
    w, V = np.linalg.eigh(K)
    return 0.5 * (V / np.sqrt(w)) @ V.T, 0.5 * (V * np.sqrt(w)) @ V.T
