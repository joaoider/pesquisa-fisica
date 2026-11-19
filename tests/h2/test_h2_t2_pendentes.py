"""H2 · T1/T2 — Próximos testes planejados (ainda sem implementação)."""
import pytest


@pytest.mark.skip(reason="Pendente: ponte H1↔H2. A entropia de emaranhamento de área, S = A/(4G), deve usar o MESMO "
                         "G induzido (Jacobson 1994; Susskind-Uglum 1994). Para escalares: S = AΛ²/(48π) ↔ 1/G = Λ²/(12π). "
                         "Vetores de gauge têm termo de contato (Kabat 1995) que é o análogo do κ = -2/3.")
def test_entropia_de_area_usa_o_mesmo_G_induzido():
    ...
