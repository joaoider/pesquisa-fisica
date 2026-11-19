"""Férmions na rede quadrada 2D num fundo de gauge U(1) com fluxo magnético (exato, redes grandes).

H = -t Σ_<ij> (e^{i A_ij} c_i^† c_j + h.c.),   fluxo por plaqueta Φ = Σ_□ A  (mod 2π).
  - Φ = 0: metal com superfície de Fermi;  Φ = π: férmions de Dirac (como a QED₃);  Φ = 2πp/q: Hofstadter.
  - Transformação de gauge U(1): c_j → e^{iθ_j} c_j, A_ij → A_ij + θ_i - θ_j. Só Φ é físico.
O campo de gauge aqui é ESTÁTICO (campo médio): é a aproximação de N grande da QED₃ / líquidos de spin U(1).
"""
from __future__ import annotations

import numpy as np


def site(x: int, y: int, L: int) -> int:
    return (x % L) + L * (y % L)


def landau_gauge_links(L: int, flux: float, twist: float = 0.0):
    """Ligações (i, j, A_ij) no calibre de Landau: A_y(x) = flux·x. Exige L·flux/2π inteiro no toro.
    twist: fase de contorno extra na direção x (evita modos zero degenerados)."""
    if not np.isclose((L * flux / (2 * np.pi)) % 1, 0) and not np.isclose((L * flux / (2 * np.pi)) % 1, 1):
        raise ValueError("fluxo incompatível com o toro: L·Φ/2π precisa ser inteiro")
    links = []
    for y in range(L):
        for x in range(L):
            i = site(x, y, L)
            links.append((i, site(x + 1, y, L), twist / L))
            links.append((i, site(x, y + 1, L), flux * x))
    return links


def hopping_matrix(L: int, links, t: float = 1.0) -> np.ndarray:
    h = np.zeros((L * L, L * L), dtype=complex)
    for i, j, A in links:
        h[i, j] += -t * np.exp(1j * A)
        h[j, i] += -t * np.exp(-1j * A)
    return h


def gauge_transform_links(links, theta: np.ndarray):
    return [(i, j, A + theta[i] - theta[j]) for i, j, A in links]


def random_flux_links(L: int, rng, strength: float = np.pi):
    """Fundo desordenado: A_ij aleatório uniforme em [-strength, strength] em cada ligação."""
    links = []
    for y in range(L):
        for x in range(L):
            i = site(x, y, L)
            links.append((i, site(x + 1, y, L), rng.uniform(-strength, strength)))
            links.append((i, site(x, y + 1, L), rng.uniform(-strength, strength)))
    return links


def plaquette_flux(L: int, links) -> np.ndarray:
    A = {}
    for i, j, a in links:
        A[(i, j)] = a
        A[(j, i)] = -a
    out = np.zeros(L * L)
    for y in range(L):
        for x in range(L):
            s0, s1, s2, s3 = site(x, y, L), site(x + 1, y, L), site(x + 1, y + 1, L), site(x, y + 1, L)
            out[s0] = np.angle(np.exp(1j * (A[(s0, s1)] + A[(s1, s2)] + A[(s2, s3)] + A[(s3, s0)])))
    return out


def ground_state_correlations(h: np.ndarray, filling: float = 0.5) -> tuple[np.ndarray, float]:
    """C_ij = <c_i^† c_j> para os N·filling modos mais baixos; devolve também o gap no nível de Fermi."""
    eps, U = np.linalg.eigh(h)
    n_occ = int(round(filling * len(eps)))
    occ = U[:, :n_occ]
    C = (occ.conj() @ occ.T)
    return C, float(eps[n_occ] - eps[n_occ - 1])


def _h2(p):
    p = np.clip(p, 1e-15, 1 - 1e-15)
    return -(p * np.log(p) + (1 - p) * np.log(1 - p))


def site_mutual_information(C: np.ndarray) -> np.ndarray:
    """I(i:j) entre sítios, para estados gaussianos de férmions complexos (vetorizado)."""
    n = len(C)
    d = np.real(np.diag(C))
    s1 = _h2(d)
    iu, ju = np.triu_indices(n, 1)
    a, b, c = d[iu], d[ju], np.abs(C[iu, ju]) ** 2
    tr, det = a + b, a * b - c
    disc = np.sqrt(np.maximum(tr**2 / 4 - det, 0))
    l1, l2 = tr / 2 + disc, tr / 2 - disc
    I = np.zeros((n, n))
    I[iu, ju] = np.maximum(s1[iu] + s1[ju] - _h2(l1) - _h2(l2), 0)
    return I + I.T


def torus_distances(L: int) -> np.ndarray:
    xy = np.array([(s % L, s // L) for s in range(L * L)])
    d = np.abs(xy[:, None, :] - xy[None, :, :])
    d = np.minimum(d, L - d)
    return np.sqrt((d**2).sum(-1))
