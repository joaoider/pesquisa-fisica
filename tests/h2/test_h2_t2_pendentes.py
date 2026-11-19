"""H2 · T1/T2 — Próximos testes planejados (ainda sem implementação)."""
import pytest


@pytest.mark.skip(reason="Pendente: ponte H1↔H2. A entropia de emaranhamento de área, S = A/(4G), deve usar o MESMO "
                         "G induzido (Jacobson 1994; Susskind-Uglum 1994). Para escalares: S = AΛ²/(48π) ↔ 1/G = Λ²/(12π). "
                         "Vetores de gauge têm termo de contato (Kabat 1995) que é o análogo do κ = -2/3.")
def test_entropia_de_area_usa_o_mesmo_G_induzido():
    ...


@pytest.mark.skip(reason="Pendente: primeira lei do emaranhamento → equações de Einstein linearizadas "
                         "(Faulkner, Guica, Hartman, Myers, Van Raamsdonk 2014), num modelo holográfico de brinquedo.")
def test_einstein_linearizado_do_emaranhamento():
    ...


@pytest.mark.skip(reason="Pendente: trocar matéria livre por campos de gauge na rede (Z2 → U(1) → SU(N)) "
                         "e repetir a reconstrução de geometria — o teste que liga H2 às '3 forças'.")
def test_geometria_a_partir_de_teoria_de_gauge_na_rede():
    ...


@pytest.mark.skip(reason="Pendente: assinatura lorentziana — de onde vem o sinal (-,+,+,+)? Combinar Page-Wootters "
                         "com o espaço emergente e testar causalidade (cone de luz de Lieb-Robinson).")
def test_assinatura_lorentziana_e_cone_de_luz():
    ...
