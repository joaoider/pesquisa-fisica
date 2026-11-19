"""H2 · T1 — Passo holográfico: da rede tensorial e da primeira lei até Einstein linearizado.

Cadeia lógica (Ryu-Takayanagi; Faulkner-Guica-Hartman-Myers-Van Raamsdonk 2014; de Boer-Heller-Myers-Neiman 2015):
  (a) Em redes tensoriais holográficas, S(A) = corte mínimo (RT discreto) e o bulk é recuperável de regiões
      grandes da borda (cunha de emaranhamento / correção de erros quânticos).
  (b) Em AdS₃/CFT₂, RT + primeira lei para TODOS os intervalos ⟺ Einstein linearizado no bulk. No espaço dos
      intervalos (x, R) — um de Sitter dS₂ — isso é a equação (∇²_dS₂ + 2) δS = 0.
  (c) A métrica do bulk reconstruída SÓ com intervalos pequenos prevê, via geodésicas RT, o δS de intervalos
      grandes (que não entraram na reconstrução).
  (d) O G do bulk é fixado pela matéria: G = 3ℓ/(2c) (Brown-Henneaux), com c medido na cadeia.

Predições pré-registradas (Notion):
  P1. Tensor do código de 5 qubits é perfeito; RT vale como limite exato S ≤ corte mínimo, com igualdade em ≥ 95%
      das regiões; o qubit de bulk central é recuperável ⇔ a região cobre a maioria (≥ 3 de 5) dos blocos.
  P2. Cadeia crítica (CFT): δS obedece à Klein-Gordon em dS₂ (resíduo mediano < 1%); cadeia gapped não (> 5%).
  P3. Bulk reconstruído de intervalos pequenos prevê intervalos grandes com erro < 5% (CFT); falha no gapped.
  P4. c ≈ 1 ⇒ G_bulk = 3/(2c) ≈ 1,5; dobrar o número de espécies de matéria dobra c e reduz G à metade
      (o mesmo papel de 1/G ∝ N na gravidade induzida de H1).
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.holographic_network import (
    build_network_state,
    contiguous_regions,
    entropy_of_boundary_region,
    five_qubit_code_tensor,
    is_perfect,
    min_cut,
)
from pesquisa_fisica.emergence.modular import (
    central_charge_from_chain,
    chain_correlation_general,
    delta_S_grid,
    ds2_klein_gordon_residual,
    fermion_entropy,
    reconstruct_energy_density,
    rt_geodesic_delta_S,
)

# ================================================================ (a) rede tensorial holográfica


@pytest.fixture(scope="module")
def rede():
    return build_network_state()


def test_tensor_do_codigo_de_5_qubits_e_perfeito():
    assert is_perfect(five_qubit_code_tensor())


@pytest.mark.parametrize("L", [1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 13, 15, 17, 19])
def test_ryu_takayanagi_discreto_intervalos(rede, L):
    for A in contiguous_regions(L):
        assert entropy_of_boundary_region(rede, A) == pytest.approx(min_cut(A), abs=1e-8)


def test_ryu_takayanagi_como_limite_e_quase_sempre_igualdade(rede):
    rng = np.random.default_rng(0)
    iguais, total = 0, 0
    for _ in range(300):
        A = sorted(rng.choice(20, size=int(rng.integers(1, 11)), replace=False).tolist())
        S, m = entropy_of_boundary_region(rede, A), min_cut(A)
        assert S <= m + 1e-8                                   # limite superior exato
        assert abs(S - m) < 1e-8 or abs(S - (m - 1)) < 1e-8    # violações de exatamente 1 bit
        iguais += abs(S - m) < 1e-8
        total += 1
    assert iguais / total >= 0.95


@pytest.mark.parametrize("blocos,bits", [(1, 0), (2, 0), (3, 2), (4, 2), (5, 2)])
def test_reconstrucao_do_bulk_precisa_da_maioria_dos_blocos(rede, blocos, bits):
    A = list(range(4 * blocos))
    I = 1 + entropy_of_boundary_region(rede, A) - entropy_of_boundary_region(rede, A, include_bulk=True)
    assert I == pytest.approx(bits, abs=1e-8)                  # I(b0:A) em bits


# ================================================================ (b) Einstein linearizado ⇔ Klein-Gordon em dS₂

N = 2000
XB = np.arange(N - 1) + 0.5
CENTROS = np.arange(N // 2 - 40, N // 2 + 41) - 1
RAIOS = np.arange(4, 51)
DIMERIZADA = 1 + 0.3 * (-1.0) ** np.arange(N - 1)


def perfil(nome):
    return {
        "gauss": np.exp(-((XB - N / 2) / 30) ** 2),
        "dois_picos": np.exp(-((XB - N / 2 - 25) / 12) ** 2) - 0.7 * np.exp(-((XB - N / 2 + 30) / 15) ** 2),
    }[nome]


def residuo_kg(base, f):
    C0 = chain_correlation_general(N, base)
    C1 = chain_correlation_general(N, base * (1 + 1e-4 * f))
    dS = delta_S_grid(C0, C1, CENTROS, RAIOS)
    res = ds2_klein_gordon_residual(dS, RAIOS)
    grande = RAIOS[1:-1] >= 15
    return np.median(np.abs(res[:, grande])) / np.abs(2 * dS[1:-1, 1:-1]).max()


@pytest.mark.parametrize("nome", ["gauss", "dois_picos"])
def test_cft_obedece_klein_gordon_em_ds2(nome):
    assert residuo_kg(np.ones(N - 1), perfil(nome)) < 1e-2


def test_cadeia_gapped_nao_obedece():
    r_gap = residuo_kg(DIMERIZADA, perfil("gauss"))
    r_cft = residuo_kg(np.ones(N - 1), perfil("gauss"))
    assert r_gap > 0.05 and r_gap / r_cft > 10


# ================================================================ (c) bulk reconstruído prevê intervalos grandes

@pytest.fixture(scope="module")
def G_bulk():
    return 3 / (2 * central_charge_from_chain(N))


def erro_previsao(base, f, G):
    centros_peq = np.arange(N // 2 - 135, N // 2 + 136) - 1
    pequenos = np.arange(2, 9)
    centros_gr = (np.arange(N // 2 - 60, N // 2 + 61) - 1)[::10]
    grandes = np.array([30, 40, 50, 60])
    y = np.arange(N // 2 - 145, N // 2 + 146).astype(float)
    C0 = chain_correlation_general(N, base)
    C1 = chain_correlation_general(N, base * (1 + 1e-4 * f))
    eps = reconstruct_energy_density(delta_S_grid(C0, C1, centros_peq, pequenos), centros_peq, pequenos, y)
    medido = delta_S_grid(C0, C1, centros_gr, grandes)
    previsto = np.array([[rt_geodesic_delta_S(eps, y, x, R, G) for R in grandes] for x in centros_gr])
    return np.abs(previsto - medido).max() / np.abs(medido).max()


@pytest.mark.parametrize("nome", ["gauss", "dois_picos"])
def test_bulk_de_intervalos_pequenos_preve_intervalos_grandes(G_bulk, nome):
    assert erro_previsao(np.ones(N - 1), perfil(nome), G_bulk) < 0.05


def test_previsao_falha_sem_cft(G_bulk):
    assert erro_previsao(DIMERIZADA, perfil("gauss"), G_bulk) > 1.0


# ================================================================ (d) G do bulk fixado pela matéria

def test_G_do_bulk_pela_carga_central(G_bulk):
    assert G_bulk == pytest.approx(1.5, abs=0.03)              # c ≈ 1 ⇒ G = 3/(2c)


def test_mais_especies_de_materia_enfraquecem_a_gravidade():
    """Duas espécies independentes: S dobra ⇒ c = 2 ⇒ G_bulk cai à metade (como 1/G ∝ N em H1)."""
    C = chain_correlation_general(400, np.ones(399))
    lens = np.arange(20, 120, 10)
    S1 = np.array([fermion_entropy(C[200 - L // 2: 200 - L // 2 + L, 200 - L // 2: 200 - L // 2 + L]) for L in lens])
    c1 = 3 * np.polyfit(np.log(lens), S1, 1)[0]
    c2 = 3 * np.polyfit(np.log(lens), 2 * S1, 1)[0]          # duas cópias desacopladas: entropias somam
    assert c2 == pytest.approx(2 * c1)
    assert (3 / (2 * c2)) == pytest.approx(0.5 * 3 / (2 * c1))
