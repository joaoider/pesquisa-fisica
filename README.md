# Pesquisa Física — gravidade emergente das três forças

**Tese:** só existem três interações fundamentais (forte, fraca e eletromagnética). A gravitação é emergente: o termo de Einstein-Hilbert e a constante G surgem das flutuações quânticas dos campos do Modelo Padrão (gravidade induzida de Sakharov, 1967; Adler-Zee).

Anotações, literatura, predições pré-registradas e diário de decisões ficam no Notion (página *Pesquisa Física — Teoria Fundamental*). Este repositório guarda os cálculos e a bateria de testes.

## Método

Cada ideia é uma hipótese que precisa sobreviver a testes cada vez mais duros:

| Nível | Pergunta | Arquivos |
|---|---|---|
| **T0** Consistência interna | A matemática fecha? | `tests/test_t0_*.py` |
| **T1** Limites conhecidos | Recupera a física estabelecida? | `tests/test_t1_*.py` |
| **T2** Vínculos experimentais | Sobrevive aos dados? | `tests/test_t2_*.py` |
| **T3** Predição nova | Diz algo novo e mensurável? | *(a criar)* |

Convenções nos testes:
- **teste normal**: resultado que a ideia precisa satisfazer;
- **`xfail(strict=True)`**: problema conhecido e aberto. Falha hoje, e isso é informação. Se passar a "passar", o pytest avisa;
- **`skip`**: teste planejado, ainda sem implementação.

A cada mudança de pressuposto, rode a suíte inteira (teste de regressão para a física).

```bash
pip install -e ".[dev]"
pytest
```

## Estrutura

```
src/pesquisa_fisica/
  geometry.py         # Christoffel, Riemann, Ricci, invariantes, □ (SymPy)
  heat_kernel.py      # coeficientes de Seeley-DeWitt a_1, a_2 + espectro exato em S^d
  induced_gravity.py  # 1/G induzido por conteúdo de campos (corte de tempo próprio)
  schemes.py          # Pauli-Villars, regularização dimensional/zeta, sensibilidade de sinal
tests/                # bateria T0–T3
notebooks/            # explorações
docs/                 # derivações e notas
```

## Resultados até agora (v0.1)

1. **Calibração (T1, passou).** O coeficiente a_1 = (1/6 − ξ)R do heat kernel foi verificado por duas rotas independentes: fórmula simbólica e soma do espectro exato do Laplaciano em S² e S⁴. Também a_2 confere. Com corte de tempo próprio Λ, um escalar real mínimo induz **1/G = Λ²/(12π)**.
2. **Sinal de G no Modelo Padrão (T0, passou, mas frágil).** Com pesos κ por campo (escalar 1/6 − ξ, Weyl +1/6, vetor −2/3, já com fantasmas), o Modelo Padrão dá Σκ = 1/6 − 4ξ, isto é, N₀ + N½ − 4N₁ = 4 + 45 − 48 = **1** em unidades de escalar.
   - Gravidade atrativa só se o acoplamento do Higgs à curvatura for **ξ < 1/24**. Higgs conforme (ξ = 1/6) dá G < 0.
   - Com neutrinos de mão direita, Σκ = 2/3 − 4ξ e ξ_c = 1/6.
   - Os bósons de gauge **sozinhos** (as "três forças" sem matéria) induzem **G < 0**. A gravidade atrativa vem da matéria, não das forças.
   - O corte necessário é Λ = √(12π) M_Pl ≈ 6,1 M_Pl.
3. **Constante cosmológica (T2, falha conhecida).** O termo a_0 induz ρ_vac ~ Λ⁴, ~10¹²⁰ acima do observado (registrado como `xfail`).

## Resultados v0.2 — dependência de esquema

4. **Pesos por spin conferidos de forma independente (T1, passou).** Somando os espectros exatos de Dirac e de 1-formas na S⁴, sem usar a fórmula de Seeley-DeWitt: tr a₁(D̸²) = −R/3 e vetor de gauge com fantasmas = −2R/3 (2 polarizações). Isso confirma os pesos escalar 1 : Weyl 1 : vetor −4. O coeficiente 1/(2π) de Λ² também bate com Visser (2002, eq. 21), a menos da convenção de sinal dele.
5. **Pauli-Villars (T0).** Para cada espécie, o sinal da contribuição é o sinal de κ, quaisquer que sejam as massas dos reguladores. Com reguladores degenerados, PV reproduz o corte de tempo próprio. **Mas o sinal total do MP depende das escalas relativas dos reguladores:** se os reguladores dos vetores forem só **~1% mais pesados** que os da matéria (razão crítica √(49/48) ≈ 1,0104; 1,041 com ν_R), G fica negativo.
6. **Regularização dimensional / zeta (T0).** Campos sem massa não induzem G nenhum. Com as massas do MP, a parte finita é |1/G| ~ 10⁴ GeV², **~34 ordens de grandeza abaixo de M_Pl²**, e depende da escala μ: G vira um contratermo, não uma predição.
7. **Veredito v0.2** (registrado como `xfail` em `test_sinal_de_G_independente_de_esquema`): o sinal de G induzido pelo MP **não é independente de esquema**. A tese, na forma "MP + 3 forças bastam", precisa de um ingrediente extra: uma física de reguladores na escala de Planck (campos pesados, cordas, estrutura discreta…) cujo espectro decide G. Essa física passa a ser parte do pressuposto.

## Próximos passos

- [x] Conferir os pesos κ e o resultado N₀ + N½ − 4N₁ (verificação espectral independente; coeficiente confere com Visser 2002)
- [x] Testar a dependência de esquema (tempo próprio, Pauli-Villars, dimensional/zeta)
- [ ] Reformular o pressuposto: que setor de reguladores (compatível com "só 3 forças") fixa G > 0 e Λ ~ M_Pl? Candidatos: torres de Kaluza-Klein (Partouche & de Vaulchier 2019), espectro de cordas, QCD-like com escala de Planck
- [ ] Conferir convenções de sinal contra Adler (1982) / Frolov-Fursaev
- [ ] Weinberg-Witten: identificar formalmente qual hipótese é violada (T0)
- [ ] Ação efetiva quadrática em h_ab: gráviton sem massa, 2 polarizações, c_g = c (T1/T2)
- [ ] T3: coeficientes de R² e Weyl² fixados pelo Modelo Padrão (a_2) e desvios de Newton em curtas distâncias

## Referências principais

- A. D. Sakharov, *Vacuum quantum fluctuations in curved space and the theory of gravitation* (1967); reimpr. Gen. Rel. Grav. 32, 365 (2000)
- M. Visser, *Sakharov's induced gravity: a modern perspective*, arXiv:gr-qc/0204062
- S. L. Adler, Rev. Mod. Phys. 54, 729 (1982)
- D. V. Vassilevich, *Heat kernel expansion: user's manual*, Phys. Rept. 388, 279 (2003)
- S. Weinberg, E. Witten, Phys. Lett. B 96, 59 (1980)
- H. Partouche, B. de Vaulchier, *Induced cosmological and Newton constants from towers of states*, arXiv:1912.09939
