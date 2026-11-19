"""H2 · T1 — SU(2) na rede em 1+1D com férmions (análogo não abeliano do Schwinger), exato.

Primeiro passo rumo às forças fraca (SU(2)) e forte (SU(3)). Partição FERMIÔNICA dos sítios
(sinais de Jordan-Wigner tratados), cada sítio = 2 modos de cor.

Predições pré-registradas (Notion):
  P1. Solver correto: J = 0 reproduz duas cópias de férmions livres; o estado fundamental é singleto de cor.
  P2. O termo elétrico de cor (Casimir S_n·S_n) é não local E não abeliano: troca de cor entre sítios
      distantes, com a mesma força a qualquer distância.
  P3. Mesmo assim, a geometria emergente do emaranhamento é a cadeia local.
  P4. Acoplamento de cor mais forte (confinamento de cor) ⇒ emaranhamento mais local, como no U(1).
"""
import numpy as np
import pytest

from pesquisa_fisica.emergence.geometry_from_entanglement import knn_graph, spearman, upper
from pesquisa_fisica.emergence.su2_lattice import SU2Chain, site_mutual_information

N = 8                                    # 8 sítios × 2 cores = 16 modos; setor: 12 870 estados
DIST = np.abs(np.subtract.outer(np.arange(N), np.arange(N))).astype(float)
CADEIA = np.isclose(DIST, 1.0)
ACOPLAMENTOS = [0.0, 0.5, 2.0, 5.0]


@pytest.fixture(scope="module")
def estados():
    out = {}
    for J in ACOPLAMENTOS:
        ch = SU2Chain(N, w=1.0, m=0.5, J=J)
        E, v = ch.ground_state()
        out[J] = (ch, E, v, site_mutual_information(ch.basis, v, 2 * N))
    return out


# ---------------------------------------------------------------- P1

def test_limite_sem_gauge_reproduz_fermions_livres(estados):
    ch, E, _, _ = estados[0.0]
    assert E == pytest.approx(ch.free_energy(), abs=1e-9)


@pytest.mark.parametrize("J", ACOPLAMENTOS)
def test_estado_fundamental_e_singleto_de_cor(estados, J):
    ch, _, v, _ = estados[J]
    assert ch.total_color_casimir(v) == pytest.approx(0.0, abs=1e-8)


# ---------------------------------------------------------------- P2

def test_termo_de_cor_troca_cor_entre_sitios_distantes():
    """S_n·S_n aplicado a |…cor0 no sítio 0…cor1 no sítio k…⟩ produz |…cor1 no 0…cor0 no k…⟩, com a
    MESMA amplitude para k = 2 e k = 5: troca de cor não abeliana que não decai com a distância."""
    ch = SU2Chain(N, J=1.0)

    def estado(occ):
        s = sum(1 << q for q in occ)
        v = np.zeros(ch.dim)
        v[ch._lookup[s]] = 1.0
        return v, s

    amps = []
    for k in (2, 5):
        outros = [2 * 7, 2 * 7 + 1, 2 * 6, 2 * 6 + 1, 2 * 4, 2 * 4 + 1] if k == 2 else \
                 [2 * 7, 2 * 7 + 1, 2 * 6, 2 * 6 + 1, 2 * 3, 2 * 3 + 1]
        v, _ = estado([0, 2 * k + 1] + outros)
        w = ch.casimir_apply(v, n=5)
        _, s_troca = estado([1, 2 * k] + outros)
        amps.append(abs(w[ch._lookup[s_troca]]))
    assert amps[0] > 0.4
    assert amps[0] == pytest.approx(amps[1])


# ---------------------------------------------------------------- P3

@pytest.mark.parametrize("J", ACOPLAMENTOS)
def test_geometria_emergente_e_a_cadeia_local(estados, J):
    I = estados[J][3]
    A = knn_graph(I, 2, mutual=True)
    assert (A & CADEIA).sum() == A.sum() == CADEIA.sum()
    assert spearman(upper(-I), upper(DIST)) > 0.9


def test_cadeia_maior_N10():
    """N = 10 (20 modos, 184 756 estados), acoplamento forte."""
    n = 10
    ch = SU2Chain(n, w=1.0, m=0.5, J=2.0)
    _, v = ch.ground_state()
    I = site_mutual_information(ch.basis, v, 2 * n)
    D = np.abs(np.subtract.outer(np.arange(n), np.arange(n))).astype(float)
    A = knn_graph(I, 2, mutual=True)
    assert (A & np.isclose(D, 1)).sum() == A.sum() == np.isclose(D, 1).sum()   # = 2(n-1): adjacência simétrica


# ---------------------------------------------------------------- P4

def test_confinamento_de_cor_deixa_o_emaranhamento_mais_local(estados):
    r = [estados[J][3][3, 4] / estados[J][3][2, 5] for J in ACOPLAMENTOS]   # vizinhos / distância 3
    assert all(a < b for a, b in zip(r, r[1:]))
    assert r[-1] / r[0] > 20
