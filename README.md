# Pesquisa Física — gravidade emergente das três forças

**Tese:** só existem três interações fundamentais (forte, fraca e eletromagnética). A gravitação é emergente: o termo de Einstein-Hilbert e a constante G surgem das flutuações quânticas dos campos do Modelo Padrão (gravidade induzida de Sakharov, 1967; Adler-Zee).

Anotações, literatura, predições pré-registradas e diário de decisões ficam no Notion (página *Pesquisa Física — Teoria Fundamental*). Este repositório guarda os cálculos e a bateria de testes.

## Método

Cada ideia é uma hipótese que precisa sobreviver a testes cada vez mais duros:

| Nível | Pergunta | Arquivos |
|---|---|---|
| **T0** Consistência interna | A matemática fecha? | `tests/test_t0_*.py`, `tests/h2/test_h2_t0_*.py` |
| **T1** Limites conhecidos | Recupera a física estabelecida? | `tests/test_t1_*.py`, `tests/h2/test_h2_t1_*.py` |
| **T2** Vínculos experimentais | Sobrevive aos dados? | `tests/test_t2_*.py` (H1; validação do código de H2 contra Martinez 2016) |
| **T3** Predição nova | Diz algo novo e mensurável? | `tests/h2/test_h2_t3_*.py` (2 tentativas, v0.26 e v0.31, ambas falharam) |

**Balanço por nível (24/11/2026):** T0 = 7 arquivos, T1 = 30, T2 = 2 (H1 + validação Schwinger/Martinez 2016, v0.30), T3 = 2 (ambos falharam: v0.26 e v0.31). O código de H2 já foi confrontado com dados de laboratório (v0.30); as afirmações específicas de H2 ainda não.

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
  emergence/          # H2: emaranhamento, geometria a partir da informação mútua, Page-Wootters
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

## H2 — o próprio espaço-tempo emerge (v0.3)

O veredito v0.2 mostrou que a versão de Sakharov (**H1**: dado um espaço-tempo, a *dinâmica* da gravidade emerge) não basta: G é fixado pelo setor ultravioleta. A hipótese original é mais forte:

> **H2:** as únicas interações fundamentais são forte, fraca e eletromagnética; **espaço e tempo** se formam a partir dessas forças e da matéria (via emaranhamento). O setor ultravioleta *é* a estrutura de onde o espaço-tempo sai.

Testes em `tests/h2/` (modelos de brinquedo de matéria livre na rede, cálculo exato de estados gaussianos):

| Teste | Resultado |
|---|---|
| T0 · Lei de área: S ∝ fronteira (2D), S = (c/3) ln L (1D) | **Passou.** c = 1,0002; S linear no perímetro com r > 0,9999 |
| T1 · Distâncias a partir da informação mútua (d = −ln I/I_max) | **Passou.** Spearman com as distâncias reais: 0,9995 (1D), 0,992 (2D), 0,994 (3D) |
| T1 · Dimensão emergente (dimensão de correlação) | **Passou.** 1,23 / 2,25 / 3,15 para redes 1D / 2D / 3D |
| T1 · Controle negativo: acoplamentos não locais aleatórios | **Passou.** Dimensão ≈ 6, sem correlação com uma grade (ρ ≈ 0) |
| T0 · Tempo emergente (Page-Wootters) | **Passou.** Universo estacionário; estado condicional obedece Schrödinger com fidelidade 1; sem emaranhamento relógio-sistema, não há tempo |

### v0.4 — geometria a partir de uma teoria de gauge: Z₂ com matéria

Modelo de gauge-Higgs de Z₂ (Fradkin-Shenker), calibre unitário, diagonalização exata numa rede 3×3 periódica (18 qubits de ligação):
H = −Δ Σ A_v − J Σ B_p − h Σ X_l − λ Σ Z_l, em que A_v é a carga de matéria, B_p o fluxo magnético, h o campo elétrico (confinamento) e λ o salto da matéria (Higgs). Arquivo: `tests/h2/test_h2_t1_gauge_z2.py`.

| Predição | Resultado |
|---|---|
| **P1.** Só a estrutura de gauge, sem flutuações (código tórico, h = λ = 0): informação mútua entre ligações exatamente zero, nenhuma geometria local | **Confirmada.** I_max < 10⁻¹²: o estado é emaranhado, mas de forma topológica (não local) |
| **P2.** Com campo elétrico e/ou matéria, os pares mais emaranhados são os vizinhos reais, em todas as fases | **Confirmada.** Precisão 1,0 (72 de 72) nas fases desconfinada, mista, confinada e de Higgs; vizinhos 25× a 2400× mais emaranhados que pares distantes |
| Perto da transição de fase a geometria "borra" | **Observado.** Contraste perto/longe cai para ~7× em (h, λ) = (0,6; 0,6) |
| **P3.** Dualidade elétrico-magnética: (h, 0) ≡ (0, h) | **Confirmada** (energia e espectro de I idênticos a 10⁻⁸) |

**Leitura física:** as "forças" sozinhas (gauge puro, sem flutuação) não geram geometria local; é preciso **flutuação de campo ou matéria** para que o emaranhamento fique local e revele a vizinhança. Isso repete, por outro caminho, o resultado de H1 (bósons de gauge sozinhos induzem G < 0) e a intuição de que "naturalmente precisaria de matéria".

**Limitações:** rede 3×3 é pequena demais para estimar a dimensão; a localidade ainda está na rede do Hamiltoniano; Z₂ é discreto (não é nenhuma das 3 forças reais).

### v0.5 — Z₂ com matéria em redes grandes: modelo de Kitaev no favo de mel

O modelo de Kitaev (2006) é exatamente uma teoria de gauge Z₂ (u_jk = ±1 nas ligações) acoplada a matéria de Majorana, com solução exata em qualquer tamanho de rede: até **L = 32 (1024 células, 2048 Majoranas)**, sem aproximação. A parte de gauge do emaranhamento é só termo de fronteira (Yao & Qi 2010). Arquivo: `tests/h2/test_h2_t1_z2_rede_grande.py`.

Novo estimador, que usa só a **ordem** da informação mútua (serve para fases gapped e críticas): grafo dos 6 parceiros mais emaranhados de cada célula, mais distâncias de grafo.

| Predição | Resultado |
|---|---|
| **P1.** Solver correto | **Passou.** Energia idêntica por duas rotas (espaço real × momentos, 10⁻¹⁰) e −0,7873 por spin (literatura: −0,7872); vórtices custam energia (Lieb) |
| **P2.** Invariância de gauge | **Passou.** I idêntica (10⁻¹⁰) após transformação de gauge aleatória |
| **P3.** Fase crítica B: o grafo emergente É a rede real | **Passou.** Precisão 1,0 em L = 14, 20, 26, 32; mesmo diâmetro (13) e mesma dimensão de crescimento de bolas; vale também com acoplamentos anisotrópicos |
| **P4.** Fase gapped A: 2D, métrica deformada | **Passou.** 6 vizinhos emergentes por célula, conectado, dimensão a ±0,3 da real; só 4/6 arestas coincidem: os dímeros fortes (z) deixam as células 50× mais emaranhadas ao longo de a₂, e o espaço emergente fica "comprimido" nessa direção |
| **P5.** Fundo de gauge desordenado | **Passou.** Precisão 0,78 e ρ = 0,89: degrada, não destrói |
| **P6.** Controle negativo | **Passou.** Grafo aleatório: diâmetro 5 (vs 13), ρ ≈ 0 |

**Leitura física:** com matéria acoplada a um campo de gauge Z₂, a geometria emergente é robusta em redes grandes. E ela é definida pelo **emaranhamento**, não pelo desenho da rede: na fase A a mesma rede dá uma métrica diferente. Isso é o comportamento esperado se o espaço é feito de emaranhamento.

### v0.6 — U(1) com férmions (eletromagnetismo na rede)

**A. Modelo de Schwinger (QED em 1+1D), U(1) dinâmico, diagonalização exata com N = 20 sítios staggered (10 células, 184 756 estados).** A lei de Gauss elimina o campo elétrico e deixa um Hamiltoniano de qubits com acoplamento coulombiano de *todos com todos*. Arquivo: `tests/h2/test_h2_t1_u1_schwinger.py`.

| Predição | Resultado |
|---|---|
| **P1.** Solver correto | **Passou.** J = 0 reproduz férmions livres (10⁻⁹); carga total nula |
| **P2.** Hamiltoniano de qubits não local | **Confirmado.** 90% dos pares acoplados; qubits a 10 sítios acoplam 9× mais que um par vizinho (potencial linear do confinamento) |
| **P3.** Geometria emergente é a cadeia local | **Passou.** Vizinhos recíprocos = cadeia exata (precisão = recall = 1) para J = 0; 0,5; 2 |
| **P4.** Gauge mais forte ⇒ emaranhamento mais local | **Passou.** Razão I(vizinhos)/I(distância 3): 42 → 289 → 9130 |

**B. Férmions em 2D num fundo U(1) com fluxo, rede 20×20 (400 sítios), exato.** O campo de gauge aqui é **estático** (campo médio da QED₃). Arquivo: `tests/h2/test_h2_t1_u1_fundo_2d.py`.

| Predição | Resultado |
|---|---|
| **P1.** Invariância de gauge U(1) | **Passou.** Mesmo fluxo e mesma I (10⁻¹⁵) após fases aleatórias e^{iθ_j} |
| **P2.** Rede emerge em todo fundo de fluxo | **Passou.** Precisão 1 para Φ = 0 (metal), π (Dirac), 2π/4 e 2π/5 (Hofstadter) |
| **P3.** Fluxo muda a escala, não a vizinhança | **Passou.** I entre vizinhos: 0,168 (Φ = 0) × 0,239 (Φ = π) |
| **P4.** Fluxo aleatório não destrói | **Passou.** Precisão 1 com desordem até ±π |
| **P5.** Controle negativo | **Passou.** Diâmetro 6 (vs 20), ρ ≈ 0 |

Observação: a meio preenchimento a informação mútua entre sítios da mesma sub-rede é exatamente zero (simetria partícula-buraco).

**Leitura física:** o resultado mais forte é o do Schwinger. O Hamiltoniano, escrito nos graus de liberdade de qubit, é não local; ainda assim, a geometria emergente do emaranhamento é a cadeia local, e o confinamento (acoplamento de gauge forte) **aumenta** essa localidade. Nesse modelo, é a força eletromagnética que "aperta" a geometria.

**Limitações:** U(1) dinâmico só em 1+1D; em 2D o fundo é estático (campo médio). QED₃ dinâmica exige DMRG/Monte Carlo.

### v0.7 — SU(2) na rede em 1+1D (análogo não abeliano do Schwinger)

Formulação de Bañuls et al. (PRX 2017): férmions staggered de 2 cores, links removidos por transformação de gauge, campo elétrico de cor dado pela carga acumulada; termo elétrico = Casimir S_n·S_n. Setor físico: meio preenchimento + singleto de cor. Diagonalização exata (Lanczos) com N = 8 (12 870 estados) e N = 10 (184 756 estados). Emaranhamento calculado com a **partição fermiônica** correta (sinais de Jordan-Wigner). Arquivo: `tests/h2/test_h2_t1_su2.py`.

| Predição | Resultado |
|---|---|
| **P1.** Solver correto | **Passou.** J = 0 = duas cópias de férmions livres (10⁻⁹); estado fundamental é singleto de cor (S² = 0) para todo J |
| **P2.** Termo de cor não local e não abeliano | **Confirmado.** O Casimir troca a cor entre sítios distantes com a mesma amplitude a distância 2 e 5 |
| **P3.** Geometria emergente é a cadeia local | **Passou.** Vizinhos recíprocos = cadeia exata para J = 0; 0,5; 2; 5 e em N = 10 |
| **P4.** Confinamento de cor ⇒ emaranhamento mais local | **Passou.** Razão I(vizinhos)/I(distância 3): 9,7 → 16,7 → 54,5 → 296 |

**Robustez do U(1):** refiz o Schwinger com a partição fermiônica. As conclusões se mantêm (cadeia exata; razão 73 → 354 → 9599), então não dependem da escolha de partição.

**Leitura física:** a força não abeliana reproduz o padrão do eletromagnetismo. O Hamiltoniano acopla e troca cor entre sítios distantes, mas o emaranhamento escolhe a vizinhança local, e o **confinamento de cor aperta a geometria**. Com U(1) e SU(2) mostrando o mesmo comportamento, ele parece ser uma propriedade geral de teorias de gauge confinantes com matéria, não um acidente de um modelo.

**Limitações:** 1+1D (a lei de Gauss elimina o campo exatamente; em 2+1D o campo de cor é dinâmico e o espaço de Hilbert cresce muito rápido); SU(2), não SU(3); redes curtas (até 10 sítios).

### v0.8 — SU(3) na rede em 1+1D (QCD 1+1D, 3 cores)

Código geral SU(N_c) (`sun_lattice.py`), com o Casimir de cor pela identidade de Fierz: S·S = ½ Σ E^{cc'}E^{c'c} − N²/(2N_c). Arquivo: `tests/h2/test_h2_t1_su3.py` (N = 6, 48 620 estados). A versão grande, **N = 8 com 2 704 156 estados**, roda no GitHub Actions: workflow manual *Cálculo pesado*, `scripts/su3_n8.py`, resultado em `results/su3_n8.json`.

| Predição | Resultado (N = 6) |
|---|---|
| **P1.** Solver correto | **Passou.** Código geral reproduz o SU(2) (10⁻¹⁴); J = 0 = 3 cópias de férmions livres; singleto de cor para todo J |
| **P2.** Cadeia emergente para todo J | **Passou** |
| **P3.** Confinamento de cor ⇒ mais local | **Passou.** Razão I(vizinhos)/I(distância 3): 53 → 73 → 212 → 2176 |

**N = 8 no GitHub Actions (2 704 156 estados; runner público 4 vCPU / 15 GB; 38 min).** [Execução](https://github.com/joaoider/pesquisa-fisica/actions/runs/44568396438), resultados em `results/su3_n8.json`, conferidos por `tests/h2/test_h2_t1_resultados_pesados.py`.

| J | Energia | Singleto (S²) | Cadeia | Spearman | I(viz.)/I(d=3) |
|---|---|---|---|---|---|
| 0 | −15,849563 (= férmions livres) | 4×10⁻²² | ✔ | 0,903 | 9,7 |
| 0,5 | −14,204158 | 5×10⁻²² | ✔ | 0,966 | 24,6 |
| 2 | −11,121905 | 1×10⁻²² | ✔ | 0,966 | 169 |
| 5 | −8,687099 | 9×10⁻²³ | ✔ | 0,966 | 2132 |

Veredito: **P1, P2 e P3 confirmadas em N = 8.** Informação mútua média por distância (1…7): em J = 0 oscila (1,47; 0,098; **0,104**; 0,022; **0,028**; 0,005; **0,018**); em J = 0,5 decai (1,35; 0,13; 0,045; 0,012; 0,004; 0,001; 0,0007); em J = 5 cai a ~0 a partir da distância 3.

**Achado (reclassificado em 19/11/2026 como reprodução, ver "Checagem de novidade"):** sem a força (J = 0), a informação mútua a distância 3 **supera** a de distância 2 (0,135 × 0,099): as duas sub-redes de férmions livres se correlacionam entre si, e a ordem por distância se inverte (Spearman 0,77). Em N = 4 o efeito chega a ligar as pontas da cadeia. **Qualquer J > 0 restaura a ordem** (Spearman ≥ 0,91). Ou seja, a matéria sozinha não define uma métrica bem-comportada, e a força de cor é o que a organiza.

### v0.9 — pendências estruturais de H2: localidade, cone de luz e o caminho para Einstein

**1. Localidade a partir do espectro** (Cotler, Penington & Ranard 2019): `tests/h2/test_h2_t0_localidade_do_espectro.py`. Usando **só a lista de autovalores**, ajusta-se um ansatz de cadeia local (12 parâmetros) a espectros de 32 e 64 níveis.

| Caso | Erro do ajuste |
|---|---|
| Espectro de Hamiltoniano local (n = 5 e 6) | **~10⁻¹⁴** (ajuste exato) |
| Espectro genérico (GUE) | ~7×10⁻³ (falha) |
| Espectro de Hamiltoniano não local | ~10⁻² a 3×10⁻² (falha) |

Onze ordens de grandeza separam "existe fatoração local" de "não existe": o espectro sabe se o sistema pode ser lido como uma cadeia de vizinhos. **Inconclusivo:** com ansatz local *sem* invariância de translação (25 parâmetros, 128 níveis), o otimizador fica preso em mínimos locais (erro ~1,5×10⁻³ mesmo para alvo local). Encontrar a fatoração pode ser difícil mesmo quando ela existe.

**2. Cone de luz, estrutura causal e isotropia de Lorentz:** `tests/h2/test_h2_t1_cone_de_luz.py`.

| Predição | Resultado |
|---|---|
| Velocidade máxima finita (Lieb-Robinson) | Férmions livres, 401 sítios: frente linear com v = 2,13 (máx. de grupo 2); fora do cone, amplitude ~10⁻⁷ |
| Estrutura causal (−,+) | Só 7×10⁻⁶ do peso em s² = (vt)² − r² < 0 |
| Vale com interação | Cadeia de Heisenberg: chegada linear, v ≈ 6,7; sítio distante intocado em t = 0,3 (< 10⁻³) |
| Isotropia de Lorentz emergente | Cone de Dirac (fluxo π): anisotropia ∝ q² (1,3% → 0,02%) e E/q → v = 2 (velocidade-limite); metal: 44% de anisotropia e E ∝ q² (sem velocidade-limite) |

**3. Rumo a Einstein, com os ingredientes de Jacobson e Faulkner et al.:** `tests/h2/test_h2_t1_primeira_lei.py` (cadeia crítica de férmions, CFT c = 1).

| Predição | Resultado |
|---|---|
| Hamiltoniano de emaranhamento local e geométrico (Bisognano-Wichmann, "temperatura de Unruh" β(x) = 2πx(ℓ−x)/(ℓv)) | **Passou.** R² = 0,998; coeficiente −1,04 a −1,05 (correção de rede de 5%); alcance par nulo; alcance 3 ≈ 2–3% |
| Primeira lei δS = δ⟨K⟩ para perturbações da "métrica" (salto modulado) | **Passou.** Razão 0,990–0,999 com K exato e com K geométrico; linear em ε |
| Reescala uniforme não muda nada | **Passou** (δS = δK = 0) |

Nota de método: perturbações de potencial local não servem para testar a primeira lei a meio preenchimento (por simetria partícula-buraco só mudam a energia de ligação em 2ª ordem); com N = 400 apareceu um desvio de tamanho finito (razão até 0,63), resolvido com N = 2000.

**Pendente:** o passo holográfico, que é impor δS = δ⟨K⟩ para *todas* as regiões e reconstruir um bulk com dimensão extra (⇒ Einstein linearizado). Exige um modelo de brinquedo holográfico (redes tensoriais HaPPY/MERA).

### v0.10 — passo holográfico: da rede tensorial até Einstein linearizado

Arquivo: `tests/h2/test_h2_t1_holografia_einstein.py`. A cadeia lógica é: RT (Ryu-Takayanagi), mais a primeira lei para **todos** os intervalos, equivale a Einstein linearizado no bulk (Faulkner et al. 2014). Em AdS₃/CFT₂ isso é a equação **(∇²_dS₂ + 2) δS = 0** no espaço dos intervalos (x, R), que é um de Sitter dS₂ (de Boer, Heller, Myers & Neiman 2015).

| Predição | Resultado |
|---|---|
| **(a)** Rede de tensores perfeitos (código de 5 qubits, árvore de 6 tensores, 20 qubits de borda + 1 de bulk) | Tensor perfeito ✔. **RT discreto exato** S(A) = corte mínimo em todos os intervalos contíguos testados; em regiões arbitrárias, S ≤ corte mínimo sempre e igualdade em 98% (violações de exatamente 1 bit) |
| **(a)** Cunha de emaranhamento / correção de erros | Qubit de bulk recuperável (I = 2 bits) ⇔ a região cobre ≥ 3 de 5 blocos; 1 ou 2 blocos: I = 0 |
| **(b)** Klein-Gordon em dS₂ para δS calculado do estado quântico (cadeia de 2000 sítios, 81 centros × 47 raios) | **CFT: resíduo mediano 0,2–0,3% (máx. ~1%)**. Gapped (dimerizada): 10–56%, ≥ 40× pior |
| **(c)** Bulk AdS₃ reconstruído **só de intervalos pequenos** (R ≤ 8) prevê, via geodésicas RT numéricas, intervalos grandes (R = 30–60) | **CFT: erro máximo 0,2–3,7%** (gaussiana, dois picos, onda). Gapped: erro de 17× (falha) |
| **(d)** G do bulk fixado pela matéria (Brown-Henneaux G = 3ℓ/2c) | c medido = 0,994 ⇒ G = 1,51; duas espécies ⇒ c = 2 ⇒ G à metade (o mesmo que 1/G ∝ N na H1) |

**Leitura:** a cadeia fecha. Matéria crítica na rede tem Hamiltoniano de emaranhamento geométrico (v0.9) e primeira lei (v0.9). Para todos os intervalos, seu δS obedece à equação que, holograficamente, é Einstein linearizado (b). Um bulk reconstruído de dados locais prevê o emaranhamento em escalas grandes (c). E a constante de Newton desse bulk é fixada pela quantidade de matéria (d). Sem criticalidade (cadeia gapped), nada disso funciona.

**Bugs encontrados e corrigidos no caminho:** fator 2 duplicado na reconstrução (1/v contado duas vezes) e cobertura insuficiente dos intervalos pequenos para perfis largos.

**Limites honestos:** AdS₃ não tem grávitons locais, então Einstein linearizado em 3D é "fraco" (a geometria fica fixada pelo tensor de energia-momento da borda). O caso não trivial é d ≥ 3 (bulk ≥ 4D), que exige CFTs em 2+1D. A rede tensorial é pequena (árvore, sem laços). E a correspondência holográfica entra como pressuposto da ligação "(b) ⟺ Einstein": o que se testa é a parte da borda.

**Limitações honestas (v0.3):** (1) o Hamiltoniano de brinquedo já é local numa rede, então o teste mostra que o *estado* codifica a geometria, não de onde vem a localidade; (2) a matéria é livre, ainda sem as 3 forças; (3) o espaço emergente é euclidiano, sem assinatura lorentziana. Cada ponto virou um teste `skip` com o próximo passo.

## Ponte H1 ↔ H2 (v0.11): S = A/(4G) usa o mesmo G induzido?

Arquivo: `tests/test_t1_ponte_h1_h2.py`. A entropia de emaranhamento de área (método das réplicas, heat kernel no cone) e o G induzido de H1 foram calculados **no mesmo esquema** (corte de tempo próprio).

| Predição | Resultado |
|---|---|
| **P1.** Termo da ponta do cone (1/12)(1/α − α) | **Passou.** Espectro exato (zeros de Bessel) em cones α = 1/2, 1/3, 1/4: erro ~10⁻⁵ |
| **P2.** Escalar mínimo: S/A = 1/(4G_ind) | **Passou, exatamente.** S/A = Λ²/48π = A/4G de H1 (10⁻⁸) |
| **P3.** Férmions: mesmo peso | **Passou.** c(Dirac) = 1,027 ≈ c(escalar) = 0,985 em 1+1D na rede; com 2 componentes transversas, Weyl = 1 unidade de escalar, igual ao peso em 1/G |
| **P4.** Escalar com ξ ≠ 0 | Diferença −6ξ (termo de Wald): o escalar conforme tem entropia mas **não induz G** |
| **P4.** Vetores de gauge | **Não fecha só com o bulk.** Contagem de estados +2 por vetor; 1/G dá −4. Falta −6 por vetor = termo de contato de Kabat = **modos de borda** do campo de gauge (Donnelly-Wall) |
| **P5.** Modelo Padrão inteiro | Sem modos de borda: **S_emaranhamento = 73 × A/(4G_ind)** (`xfail`). Com modos de borda: **razão = 1 exatamente** (com ou sem ν_R, ξ pequeno) |

**Leitura:** a ponte fecha, mas só com uma peça específica. As partes de matéria (escalares mínimos e férmions) têm entropia de emaranhamento exatamente igual a A/(4G_ind). Os campos das forças não: o que fecha a conta são os **modos de borda dos campos de gauge** na superfície de emaranhamento, com contribuição negativa (−6 por vetor). Para a tese, isso significa que, se a gravidade emerge das 3 forças, a entropia de Bekenstein-Hawking não é só "contagem de estados de matéria": ela inclui os graus de liberdade de borda das forças. Esse é o mesmo sinal negativo dos vetores em G (H1) e o mesmo papel da estrutura de gauge que não gera geometria sozinha (código tórico, v0.4).

**v0.12: o termo de contato agora é calculado, não citado.** No mesmo esquema da ponte (heat kernel + réplicas), calculei o espectro exato de uma **1-forma** no disco cônico. Ao dar a volta na ponta, o vetor gira de 2πα, então as componentes V_x ± iV_y carregam índice angular ν = m/α ± 1 (zeros de Bessel exatos para α = 1/n). Decompondo a constante do heat kernel em borda (∝ α) e ponta [(1/α − α) e (1 − α)]:

| Quantidade (unidades de escalar) | Calculado |
|---|---|
| Escalar: coeficiente de holonomia q | 2×10⁻⁵ (≈ 0) |
| 1-forma: coeficiente de holonomia q | **−0,99999** (termo −(1 − α) na ponta) |
| 1-forma no plano do cone | −4,0000 (dois escalares dariam +2) |
| Campo de gauge (2 transversas + 1-forma − 2 fantasmas) | **−4,0000** = peso do vetor em 1/G (H1) |
| Contagem física (2 polarizações) | +2,0001 |
| **Termo de contato** | **−6,0001** (Kabat) |
| Modelo Padrão com pesos calculados: S/(A/4G) | **1,000** |

**Parte universal, na rede e independente de esquema:** termos de área dependem do regulador, mas o termo logarítmico da esfera não. Para Maxwell, o resultado físico (sem modos de borda) difere do da anomalia por −1/3, que é a contribuição dos modos de borda: −2 × (coeficiente da onda parcial l = 0 de um escalar). Calculei esse coeficiente numa rede radial de Srednicki com 3000 sítios: **0,1646 a 0,1652 (1/6 com 1,3% de erro)**, o que dá −1/3 para os modos de borda.

**Limite honesto:** a parte universal identifica a contribuição dos modos de borda pelo argumento de Casini-Huerta e Donnelly-Wall (Maxwell físico = 2 escalares sem l = 0); o que calculo é o bloco l = 0 e o termo de contato no esquema da ponte. Não implementei a teoria de Maxwell completa na rede cúbica com centro elétrico.

## Checagem de novidade (19/11/2026)

Os dois achados que pareciam originais foram confrontados com a literatura e **reclassificados como corolários conhecidos**:

1. **"O confinamento aperta a geometria"** (v0.6–v0.8): explicado por *gap ⇒ correlações exponenciais ⇒ informação mútua local* (Hastings-Koma; Wolf, Verstraete, Hastings & Cirac, PRL 100, 070502, 2008). Aumentar o acoplamento de gauge aumenta o gap. Medidas de emaranhamento não distinguem teoria confinante de teoria apenas com gap (arXiv:2010.09392; arXiv:2401.01930). Leitura correta: **o gap gerado pela força** aperta a geometria. Teste discriminante planejado: comparar gap de gauge e gap sem gauge, com o mesmo gap.
2. **"Sem força, a métrica não é monotônica"** (v0.8): são as oscilações de paridade universais das entropias em cadeias críticas (Calabrese, Campostrini, Essler & Nienhuis, PRL 104, 095701, 2010), suprimidas pelo gap.

Vizinhança da tese: Szangolies (arXiv:2512.17328) é o trabalho conceitual mais próximo (grupo do MP a partir de emaranhamento de qubits; espaço-tempo da parte de lei de área; sem testes numéricos). Wang & Braunstein (arXiv:2606.27402) mostram que gauge SU(2) puro não gera geometria em estados típicos, o que é coerente com o código tórico (v0.4).

O que pode ter valor original: (i) geometria emergente de **estados fundamentais de teorias de gauge com matéria**, com controles negativos sistemáticos (precisa de busca mais funda); (ii) o teste discriminante gap × confinamento; (iii) a síntese reprodutível H1 + H2 + ponte, que tem valor metodológico e não é descoberta.

## Teste discriminante: gap × confinamento (v0.13)

Arquivos: `src/pesquisa_fisica/emergence/discriminant.py`, `tests/h2/test_h2_t1_discriminante_gap_confinamento.py`. **Hipótese nula (literatura):** com o *mesmo gap*, a geometria emergente da informação mútua é a mesma, venha o gap da força de gauge ou de uma massa sem gauge.

**U(1), Schwinger com massa de quark zero (todo o gap vem da força) × cadeia livre com massa ajustada para o mesmo gap** (N = 20, partição fermiônica). Razão I_gauge(d)/I_livre(d):

| J | gap | d = 1 | d = 2 | d = 3 | d = 4 | d = 5 | d = 6 |
|---|---|---|---|---|---|---|---|
| 0,5 | 1,053 | 1,07 | 1,26 | 1,17 | 0,96 | 0,72 | 0,51 |
| 1 | 1,591 | 1,13 | 1,17 | 0,89 | 0,58 | 0,35 | 0,20 |
| 2 | 2,543 | 1,15 | 0,93 | 0,51 | 0,24 | 0,11 | 0,04 |

**H0 rejeitada.** No mesmo gap, a força dá emaranhamento **mais local**: mais informação mútua entre vizinhos e decaimento mais rápido. O excesso cresce com o acoplamento. Taxa de decaimento ÷ gap: 1,34–1,40 com gauge contra 0,97–1,17 com massa livre.

**Controle mais duro:** igualando o comprimento de correlação (a taxa de decaimento) em vez do gap, os perfis ainda diferem. A força dá 1,1–1,5× mais informação mútua em d = 1 e 2–3× em d ≥ 2. Não é "só ξ".

**SU(2)** (quark com massa 0,5, N = 8, singleto de cor): mesmo padrão. A razão em d = 5 vai de 0,81 a 0,021 de J = 0,5 a J = 4. **E confinamento ≠ gap:** com quarks sem massa, o gap do SU(2) *diminui* com o acoplamento (méson leve, à la 't Hooft).

**Mecanismo provável (bosonização):** na teoria confinada, as excitações com gap são **mésons** (compostos bosônicos). Na livre, o portador mais leve é um férmion com metade do gap do par. Como I ~ (correlação)², no contínuo o decaimento confinado deveria ser **2×** mais rápido no mesmo gap. Na rede finita vejo cerca de 1,4×.

**Novidade (avaliação honesta):** o mecanismo é conhecido (bosonização; Florio, arXiv:2312.05298, sobre negatividade fermiônica × espectro bosônico no Schwinger). Não encontrei o teste de *gap pareado* para a geometria emergente, e o resultado contradiz, neste cenário de rede em 1+1D, a afirmação holográfica de que o emaranhamento não distingue confinamento de gap (arXiv:2010.09392, que trata de teorias holográficas de N grande, um regime diferente). Para a tese, a leitura é que **a força não é só uma fonte de gap: ela muda a natureza dos portadores de emaranhamento** (de férmions para mésons), e isso aperta a geometria.

**Limites:** redes pequenas (10 células no U(1), 8 sítios no SU(2)); só 1+1D; taxas medidas a d ≤ 6. **Próximo:** cadeias longas (DMRG) para medir a razão assintótica das taxas e testar a predição de 2 do contínuo.

## Cadeia longa: razão assintótica das taxas (v0.14) — previsão de 2 FALSIFICADA

**Pergunta.** Comparando Schwinger sem massa e férmion livre com o **mesmo gap**, a informação mútua (IM) do gauge decai mais rápido. A previsão pré-registrada (bosonização: mésons × férmions) era que a razão das taxas assintóticas tende a **2 ± 0,2** no contínuo.

**Método.** O Schwinger (Kogut-Susskind, Gauss eliminado) foi resolvido por DMRG (TeNPy, MPO de dimensão 5, χ ≤ 140) até **N = 96**, em `emergence/schwinger_dmrg.py`. O férmion livre staggered foi resolvido exatamente (gaussiano), com a massa casada pelo gap. Foram medidas três quantidades:
- a correlação **carregada** ⟨c†ᵢ cⱼ⟩, com string de Jordan-Wigner;
- a correlação de **densidade** ⟨QᵢQⱼ⟩_c;
- a IM entre dois sítios, nas partições fermiônica e de spin.

As taxas vêm de um ajuste `A d^(−α) e^(−κd)` (o ajuste exponencial puro superestimava a taxa livre em até 30%). A referência livre usa a taxa exata κ = 2·asinh(m/2).

**Validação.** O MPO reproduz a diagonalização exata (energia e gap, com erro abaixo de 1e-7). Com J = 0, o DMRG coincide com a solução gaussiana até 1e-11.

| J | gap (N=96) | κ_IM gauge / gap | κ_IM livre / gap | **razão** |
|---|---|---|---|---|
| 1,0 | 1,570 | 0,716 | 0,488 | **1,46** |
| 0,5 | 1,017 | 0,736 | 0,494 | **1,48** |
| 0,25 | 0,676 | 0,750 | 0,495 | **1,50** |
| 0,1 | 0,406 | 0,771 | 0,493 | **1,53** |

A razão é a mesma em N = 64 e 96 (diferença < 0,03). Extrapolada para J → 0, fica em ≈ **1,55–1,6**, fora de 2 ± 0,2.

| Previsão | Resultado |
|---|---|
| Razão gauge/livre → 2 no contínuo | **Falhou.** Converge a ≈ 1,5–1,6 |
| P1: no gauge, a correlação carregada decai mais rápido que a de densidade (tensão de corda) | **Falhou.** κ_D/κ_C ≈ 1,4: a carregada decai *mais devagar* |
| P2: no livre, κ_D = 2κ_C (Wick) | **Passou.** 2,00 |

**Leitura.**
- A premissa errada era que a IM do gauge seria dominada pela densidade de mésons. Ela é dominada pela correlação carregada ao quadrado, como no caso livre, mas essa correlação decai a ≈ 0,75·gap/v em vez de 0,5·gap/v.
- Em **m = 0 o Schwinger blinda cargas** (não confina): não há tensão de corda acelerando o decaimento carregado. Isso explica a falha de P1.
- O efeito robusto continua de pé: o gauge torna o emaranhamento ~1,5× mais local que um férmion livre de mesmo gap. O fator, porém, não é o 2 da bosonização ingênua.
- **Próximo teste natural:** m ≠ 0, onde há confinamento linear (σ ∝ m·g). Lá, P1 deveria passar.

## Regime confinante (v0.15): três previsões falsificadas e uma lei de escala

**Pergunta.** Na v0.14, com m = 0, o modelo gauge deixa o emaranhamento ~1,5× mais local que o férmion livre de mesmo gap. Em m = 0 o Schwinger *blinda* as cargas. Com m ≠ 0 surge a tensão de corda: o fluxo elétrico custa ≈ J por elo. Previsões pré-registradas no Notion antes de rodar:
- **P1:** a taxa local κ_C(d) da correlação carregada cresce com a distância (corda);
- **P2:** essa subida satura em d* ≈ gap/J (quebra de corda);
- **P3:** a razão de IM gauge/livre sobe acima de 1,5.

**Método.** DMRG com N = 64 e χ ≤ 160 (precisão SVD de 1e-14), com o mesmo casamento de gap da v0.14. Os pontos rodados foram (m, J) = (0; 0,25; 0,5; 1 | 0,5), (0,5 | 1; 2). O diagnóstico exploratório com quarks pesados (m = 1, J = 0,1 e 0,25; corda de até ~21 sítios) foi feito *depois* e não estava pré-registrado.

| m | J | x = m/√J | razão IM gauge/livre | κ_C(d) |
|---|---|---|---|---|
| 0 | 0,5 | 0 | **1,49** | constante |
| 0,25 | 0,5 | 0,35 | **1,36** | constante |
| 0,5 | 2,0 | 0,35 | **1,35** | constante |
| 0,5 | 1,0 | 0,50 | **1,33** | constante |
| 0,5 | 0,5 | 0,71 | **1,28** | constante |
| 1,0 | 0,5 | 1,41 | **1,20** | constante |
| 1,0 | 0,25 | 2,0 | **1,14** | cai ~7% e satura |
| 1,0 | 0,1 | 3,2 | **1,07** | cai ~10% e satura |

| Previsão | Resultado |
|---|---|
| P1: κ_C(d) cresce (corda) | **Falhou.** A maior subida é < 3% em todos os pontos com m ≥ 0,5. Com quarks pesados, κ_C(d) *cai* |
| P2: saturação em d* ≈ gap/J | **Falhou.** Não há subida para saturar |
| P3: razão > 1,5 com m > 0 | **Falhou, com sinal oposto.** A razão *cai* de 1,49 para 1,07 |

**Achado (não previsto).**
- A razão depende de um único parâmetro, x = m/√J. No contínuo, J = g²a/2 e w = 1/(2a), então x ∝ m/g.
- Pontos com J diferentes e o mesmo x colapsam: 1,36 contra 1,35.
- A razão decresce monotonicamente de ≈ 1,5 (x = 0, interação puramente gauge) até 1 (quarks pesados, x → ∞, regime quase livre).

**Leitura.**
- A localidade extra do emaranhamento vem do **acoplamento relativo à massa**, não do confinamento.
- O correlador estático do vácuo **não enxerga a corda**: com matéria dinâmica, flutuações de pares quebram a corda virtual logo de cara (lei de perímetro), e o decaimento longo fica com os mésons.
- A confinamento aparece na **dinâmica**: supressão do cone de luz após um quench (Kormos et al., Nat. Phys. 2017), ou laços de Wilson com cargas externas estáticas. É esse o teste certo para a corda.

## Confinamento na dinâmica (v0.16): quench global e o cone de luz do espaço emergente

**Pergunta.** A v0.15 mostrou que o vácuo estático não enxerga a corda. Ela aparece quando o sistema é tirado do equilíbrio? E o efeito separa confinamento (m > 0) de blindagem (m = 0)?

**Protocolo.**
- Estado inicial: vácuo nu staggered (sem pares, sem fluxo).
- Evolução: TDVP de dois sítios (TeNPy) com H(m, J), N = 24, t ≤ 6, χ = 96.
- Referência: férmion livre exato com o **mesmo gap**.
- Observáveis: entropia de meia cadeia S(t) e IM de dois sítios I(c, c+d; t).
- Validação: o TDVP reproduz a dinâmica livre exata (erro < 1e-3). χ = 96 e χ = 160 dão razões iguais (diferença < 0,01).

| m | J | S_gauge/S_livre (t = 6) | alcance da frente de IM (sítios; livre ≥ 16) |
|---|---|---|---|
| 0 | 0,5 | **0,98** | 10 |
| 0 | 1 | **0,87** | 9 |
| 0,5 | 0,5 | **0,83** | 12 |
| 0,5 | 1 | **0,66** | 9 |
| 0,5 | 2 | **0,43** | 7 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: m > 0 ⇒ razão < 0,5 e caindo com J | **Falhou como enunciado.** Cai com J (0,83 → 0,66 → 0,43), mas só J = 2 fica abaixo de 0,5 |
| P2: m = 0 (blindagem) ⇒ razão > 0,7 | **Passou.** 0,98 e 0,87. Com J fixo, ligar m reduz a razão (0,87 → 0,66; 0,98 → 0,83) |
| P3: frente de IM do gauge com velocidade < 0,6 da livre | **Falhou.** A frente não tem velocidade constante: arranca e **trava**. Trava também em m = 0, com o mesmo alcance de m = 0,5 em J = 1 |

**Leitura.**
- A corda aparece na dinâmica, como esperado da literatura (Kormos et al. 2017; Magnifico et al. 2020): a J fixo, a massa freia o crescimento do emaranhamento em larga escala (S de meia cadeia). Isso separa confinamento de blindagem, coisa que o vácuo estático não fazia.
- A trava da **IM de dois sítios** é outra coisa. Ela depende de J e não de m, então vem da interação gauge em si, e não da corda.
- No espaço emergente, o gauge "fecha" o cone de luz da geometria local (pares de sítios) mesmo sem confinamento. A corda age sobre o emaranhamento de longa escala.

## Correlações multipartidas após o quench (v0.17): a informação "sumida" não está em blocos de 2

**Pergunta.** Na v0.16, a IM entre pares de sítios do gauge trava após o quench, mas a entropia de meia cadeia quase não muda. Hipótese: a informação está em correlações de vários sítios, como mésons de 2 sítios.

**Método.**
- Mesmo quench da v0.16 (N = 24, t = 6), agora com **gauge e livre por TDVP** na mesma partição (spin/qubits).
- Matrizes densidade de blocos disjuntos calculadas por contração direta do MPS (`block_pair_rdms`, validada contra os correlatores com erro de 1e-15).
- Medidas: IM entre blocos de 2 sítios, I(A₂:B₂), e o excesso multipartido f = [I(A₂:B₂) − Σ I(a:b)] / I(A₂:B₂).
- Os blocos são **alinhados** aos pares do vácuo ([2k, 2k+1]) ou desalinhados.

| caso | alcance sítio → bloco (gauge) | alcance sítio → bloco (livre) | f gauge | f livre |
|---|---|---|---|---|
| m = 0, J = 1 | 5 → 8 | 15 → 18 | 0,57 | 0,63 |
| m = 0,5, J = 1 | 6 → 8 | 11 → 14 | 0,75 | 0,63 |
| m = 0,5, J = 2 | 4 → 8 | 11 → 14 | 0,86 | 0,31 |

| Previsão | Resultado |
|---|---|
| P1: blocos ganham mais alcance no gauge que no livre | **Falhou** em 2 de 3 casos. Os blocos também travam em ~8 sítios |
| P2: f_gauge > f_livre nos 3 casos | **Falhou** em m = 0 (0,57 < 0,63). Passou com m > 0 |
| P3: f maior com blocos alinhados (méson) | **Falhou** em J = 2 (0,862 contra 0,870). As diferenças são pequenas em todos os casos |

**Achados (não previstos).**
1. **Não há sinal balístico no gauge.** No livre, a IM de blocos tem picos na borda do cone de luz: são os pares de quasipartículas emaranhadas do quadro de Calabrese-Cardy, por exemplo d = 12 em m = 0. No gauge, a IM decai monotonicamente, com ou sem confinamento.
2. **A fração multipartida cresce com o confinamento:** 0,57 → 0,75 → 0,86.

**Leitura.**
- A entropia de meia cadeia do gauge (S = 2,4 em m = 0) não vem de pares localizados de 2 sítios a longa distância.
- Ou ela está em blocos maiores, ou é uma forma de emaranhamento sem portador local, que nenhum par de regiões pequenas enxerga.
- O próximo passo é medir a escala diretamente: a entropia de blocos contíguos S(ℓ) em função do tamanho ℓ, e a informação tripartida I₃.

## Escala do emaranhamento após o quench (v0.18): confinamento ⇒ lei de área; blindagem mantém pares longos

**Pergunta.** A v0.17 mostrou que a entropia de meia cadeia do gauge não está em pares de regiões pequenas. Em que escala ela está?

**Método.** Mesmo quench (N = 24, t = 6, TDVP), com gauge e livre de mesmo gap na partição de spin. Medidas:
- a entropia S(ℓ) de um bloco central contíguo, para ℓ = 1..10 (SVD na forma mista, validada contra o TeNPy);
- a IM entre as regiões **externas** ao bloco, I(L:R) = S_corte(c) + S_corte(c+ℓ) − S(bloco), que mede os pares que atravessam o bloco inteiro;
- a informação tripartida I₃ de três blocos adjacentes de 2 sítios.

| caso | ℓ* (90% da saturação), gauge / livre | I(L:R) em ℓ = 8, gauge / livre | I₃, gauge / livre |
|---|---|---|---|
| m = 0, J = 1 (blindagem) | 9 / 9 | 0,79 / 1,95 (**0,40**) | −0,107 / −0,161 |
| m = 0,5, J = 1 | 7 / 8 | 0,10 / 0,89 (**0,12**) | −0,185 / −0,029 |
| m = 0,5, J = 2 | 5 / 9 | 0,004 / 0,34 (**0,01**) | −0,096 / −0,082 |

| Previsão (nos 3 casos) | Resultado |
|---|---|
| P1: S(ℓ) do gauge satura antes | **Falhou em m = 0** (empate em 9). Passou com m > 0 |
| P2: I(L:R)_gauge < 0,3 · livre | **Falhou em m = 0** (0,40). Passou com m > 0 (0,12; 0,01) |
| P3: I₃ do gauge negativa e abaixo da do livre | **Falhou em m = 0**. Passou com m > 0 |

**Leitura.** As três previsões falham exatamente no caso blindado e passam nos dois casos confinados.
- **Confinamento (m > 0)** produz emaranhamento de curto alcance e denso. Com J = 2, S(ℓ) do gauge para em ~1,05 a partir de ℓ ≈ 5 (lei de área), enquanto o livre segue crescendo até 2,4. Os pares que atravessam 8 sítios praticamente somem (1% do livre), e a informação fica multipartida (I₃ mais negativa).
- **Blindagem (m = 0)**: o Schwinger sem massa é um bóson livre massivo. Os mésons propagam balisticamente e **carregam pares longos** (I(L:R) em ℓ = 8 é 40% do livre). Esses pares são invisíveis à IM de sítios e de blocos de 2 da v0.16–v0.17, porque o méson é um objeto estendido nas variáveis de rede.
- Isso resolve o quebra-cabeça da v0.17: o emaranhamento "sumido" em m = 0 está em escala grande, carregado por objetos compostos.

**Para H2.** O espaço emergente de um sistema confinante é estritamente local, com emaranhamento de lei de área após o quench. O de um sistema blindado mantém conectividade de longo alcance, mas só em escalas maiores que o tamanho do méson. A escala em que a geometria "se forma" é o tamanho das partículas compostas.

## Coarse-graining em escala de méson (v0.19): o cone de luz só é causal nos observáveis invariantes de gauge

**Pergunta.** Na v0.18, a IM entre as regiões externas mostrou pares longos em m = 0. Eles aparecem como um cone de luz balístico na IM entre blocos do tamanho de um méson?

**Método.**
- **Evolução exata** com N = 20 (setor de carga zero, 184.756 estados), sem truncamento. Validada contra a solução gaussiana (erro 1e-13).
- Quench do vácuo nu, t ≤ 6. IM entre dois blocos de tamanho ℓ_b = 1..4, simétricos em torno do centro e separados por r = 2..8.
- Diagnóstico extra: a IM **só de cargas** (diagonal de ρ na base σᶻ), que contém apenas observáveis locais invariantes de gauge.

| Previsão pré-registrada | Resultado |
|---|---|
| P1: m = 0, ℓ_b = 4 ⇒ pico balístico em todo r, com t_pico crescente | **Falhou.** Os picos são **simultâneos** (t_pico ≈ 1,6 para r = 2..8). No livre, t_pico = 2,8 → 4,0 → 4,7 → 5,05 |
| P2: v_gauge/v_livre ∈ [0,5; 1] | **Falhou.** A frente total do gauge corre a v ≈ 3,8 (acima de 2w = 2), com chegada sublinear em r |
| P3: m = 0,5, J = 2 ⇒ I_max(r = 8) < 0,1 · livre | **Passou** (0,022) |

**Diagnóstico: a frente "rápida" é artefato de corda.**

| tempo de chegada (r = 2, 4, 6, 8) | IM total | só cargas (invariante de gauge) |
|---|---|---|
| gauge m = 0 | 0,23, 0,58, 0,84, 1,02 | **0,40, 0,93, 1,44, 2,03** |
| livre de mesmo gap | 0,23, 0,63, 1,12, 1,71 | **0,40, 0,96, 1,47, 2,12** |
| gauge m = 0,5, J = 2 | 0,23, 0,52, 0,71, 0,85 | **0,41, 1,00, 2,48, 4,90** |
| livre de mesmo gap | 0,23, 0,66, 1,17, 2,37 | 0,41, 1,02, 2,44, 3,61 |

- A IM total do gauge a longa distância é quase toda **coerente**: a parte de cargas é < 5% no pico.
- A IM total vem de operadores como σ⁺_a σ⁻_b. Com a lei de Gauss eliminada, eles carregam implicitamente a **corda de fluxo** entre os blocos. São operadores estendidos: o valor esperado deles evolui pela dinâmica local ao longo de toda a corda, sem sinal viajando. Daí os picos simultâneos e a frente acima de 2w.
- Nos **observáveis locais invariantes de gauge** (cargas), o cone de luz é causal. Em m = 0 ele coincide com o do férmion livre (v ≈ 2). Com confinamento, é freado (r = 8 chega em t = 4,9 contra 3,6 no livre), como em Kormos et al.

**Consequência metodológica para H2.** Toda a IM de teorias de gauge usada desde a v0.6 (U(1), SU(2), SU(3), v0.13–v0.18) usou partições **só de matéria com a lei de Gauss eliminada**. Essa partição mistura correlações locais com coerências de corda. É a mesma questão de "qual álgebra?" que apareceu no termo de contato de Kabat (v0.12): o emaranhamento de campos de gauge depende da escolha do centro elétrico. Os resultados qualitativos de confinamento sobrevivem (P3 e a frente freada nas cargas), mas os números de localidade precisam ser refeitos numa partição invariante de gauge.

## Localidade em álgebras invariantes de gauge (v0.20): a causalidade volta; a conclusão qualitativa de H2 sobrevive

**Pergunta.** A v0.19 mostrou que a IM na partição de qubits (lei de Gauss eliminada) mistura coerências de corda. O que acontece com uma álgebra invariante de gauge?

**Álgebras comparadas** (Casini, Huerta & Rosabal 2014), em `emergence/gauge_algebra.py`:
- **qubit:** a álgebra completa da região (a usada até a v0.19). Inclui implicitamente linhas de Wilson entre regiões.
- **inv:** operadores que comutam com a carga de cada região. Não há linhas de Wilson entre regiões; ρ fica decomposto por setores Q_A, Q_B.
- **ce (centro elétrico):** inv mais os fluxos de fronteira. Em 1D aberto, o fluxo à esquerda de R é a carga total à esquerda (variável clássica).
- **carga:** só a diagonal na base de carga.

**Quench exato** (N = 20, m = 0, J = 1, blocos de 4 sítios):

| álgebra | t_pico (r = 2, 4, 6, 8) | velocidade da frente |
|---|---|---|
| qubit | 1,61 · 1,60 · 1,58 · 1,59 (simultâneos) | 3,77 |
| **inv** | **2,46 · 4,54 · 5,77 · > 6** | **1,89** |
| ce | 1,92 · 2,08 · 5,63 · > 6 | 1,76 |
| carga | 1,40 · 1,91 · 5,66 · 5,60 | 1,85 |

**Vácuo** (N = 16, m = 0,5; regiões = células de 2 sítios), razão I(d = 1)/I(d = 3):

| J | qubit | **inv** | ce | carga |
|---|---|---|---|---|
| 0 | 42 | **2,4·10³** | 6,1·10² | 2,6·10³ |
| 0,5 | 289 | **3,0·10⁴** | 1,8·10⁴ | 3,3·10⁴ |
| 2 | 9,1·10³ | **7,4·10⁶** | 5,4·10⁶ | 7,7·10⁶ |
| 5 | 7,4·10⁵ | 1,5·10¹⁰ * | 1,0·10¹⁰ * | 1,5·10¹⁰ * |

\* Em J = 5, I(d = 3) ~ 10⁻¹¹ está no piso numérico.

| Previsão pré-registrada | Resultado |
|---|---|
| P1: na álgebra inv a frente é causal (t_pico crescente, v ∈ [1,5; 2,3]) | **Passou** (v = 1,89) |
| P2: o centro elétrico é mais rápido (v_ce > v_inv + 0,5) | **Falhou** (v_ce = 1,76) |
| P3: no vácuo, a razão de localidade cresce com J na álgebra inv | **Passou** (0 → 0,5 → 2) |

**Leitura.**
- Remover as linhas de Wilson entre regiões **restaura a causalidade**: a frente volta a v ≈ 2, igual à do férmion livre.
- O centro elétrico **também é causal**. Meu raciocínio de que o fluxo de fronteira levaria informação de Coulomb instantânea estava errado: o fluxo só muda por correntes locais que atravessam a fronteira.
- A conclusão qualitativa de H2 (v0.6–v0.8: acoplar o gauge deixa o emaranhamento mais local) **sobrevive** na álgebra invariante.
- Os **números mudam por ordens de grandeza**. Em J = 0 a razão é 42 em qubits e 2.400 em inv. Boa parte da "não localidade" de base vinha de coerências de hopping/corda, e não de correlações invariantes.
- Na álgebra invariante, quase toda a IM é de cargas: inv ≈ carga. A parte quântica dentro das células é pequena.

**Para H2:** a geometria emergente de uma teoria de gauge precisa ser definida numa álgebra invariante. Com ela, o espaço emergente é causal e local, e mais local com acoplamento forte. É a mesma lição do termo de contato de Kabat (v0.12): em teoria de gauge, "emaranhamento" só tem sentido junto com a escolha da álgebra.

## Razão gauge/livre na álgebra invariante (v0.21): o fator 1,5 era artefato; a mesmo gap, gauge = livre

**Pergunta.** A v0.14 mediu que o Schwinger deixa o emaranhamento ≈ 1,5× mais local que o férmion livre de mesmo gap, e a v0.15 transformou isso numa lei de escala em m/g. As duas medidas foram feitas em partições de matéria (qubits/fermiônica). O efeito sobrevive na álgebra invariante (v0.20)?

**Argumento antes de rodar.** Na álgebra invariante só restam correlações **neutras**. Nos dois modelos, elas decaem com a excitação neutra mais leve, cuja massa é justamente o gap que está sendo casado. Previsão: razão ≈ 1.

**Método.** DMRG com N = 64 para o Schwinger (m = 0) e para o livre de mesmo gap, na mesma partição. As regiões são células de 2 sítios, e a IM entre células foi calculada nas álgebras qubit, inv e carga (validada contra a diagonalização exata, erro 1e-8). Taxa por sítio vinda do ajuste A·d^(−α)·e^(−κd).

| J | gap | **razão inv** | κ_inv/gap (gauge; livre) | razão carga | razão qubit (controle) |
|---|---|---|---|---|---|
| 0,25 | 0,680 | **1,03** | 0,97; 0,94 | 1,02 | 1,48 |
| 0,5 | 1,019 | **1,04** | 0,99; 0,95 | 1,04 | 1,48 |
| 1,0 | 1,572 | **1,07** | 1,00; 0,94 | 1,06 | 1,46 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: razão inv ∈ [0,85; 1,15] | **Passou** (1,03–1,07) |
| P2: κ/gap ∈ [0,8; 1,2] nos dois modelos | **Passou** (0,94–1,00) |
| P3: controle de qubits ∈ [1,3; 1,7] (reproduz a v0.14) | **Passou** (1,46–1,48) |

**Leitura.**
- O fator ≈ 1,5 da v0.14 e a lei de escala em m/g da v0.15 vinham de **correlações carregadas**, isto é, de linhas de Wilson implícitas.
- Na álgebra invariante, **a mesmo gap o gauge é indistinguível do livre**. Há um resíduo de 3–7%, que cresce levemente com J (efeito de rede ainda não isolado).
- Toda a localização do emaranhamento pelo gauge vem da **geração de gap** (a massa de Schwinger μ = g/√π), e não de algum efeito adicional da força.

**Consequência para H2.**
- A versão forte ("o gauge deixa o espaço emergente mais local do que a massa sozinha explicaria") está **refutada** na álgebra invariante.
- Sobrevive a versão fraca: o acoplamento gauge gera massa, e a massa localiza o emaranhamento. Isso é genérico, e qualquer teoria com gap faz o mesmo.
- O papel específico do gauge no espaço emergente tem de ser procurado em outro lugar: no termo de contato e no centro elétrico (v0.11–v0.12), na dinâmica confinante (v0.16) ou na estrutura de lei de área (v0.18).

## Confinamento na dinâmica na álgebra invariante (v0.22): o efeito é real, e o gauge deixa uma assinatura multipartida

**Pergunta.** Depois da v0.21 (no vácuo, a mesmo gap, gauge = livre na álgebra invariante), os efeitos dinâmicos do confinamento das v0.16–v0.18 também eram artefato de partição?

**Teorema (P0).** O estado tem carga total definida. Então, para **uma região só**, ρ_R comuta com Q_R, e S_inv(R) = S_qubit(R). Verificado numericamente com diferença < 1e-14. Consequência: a supressão de S(t) pelo confinamento (v0.16) e a lei de área de S(ℓ) (v0.18) **sobrevivem automaticamente**. Só as quantidades com regiões disjuntas precisavam ser refeitas.

**Método.** Quench exato (N = 20), t ≤ 6, com gauge e livre de mesmo gap. Medidas:
- I(L:R) entre as regiões externas a um bloco central de 8 sítios;
- I₃ de três blocos adjacentes de 2 sítios.

As álgebras conjuntas são geradas pelas regiões, sem linhas de Wilson entre elas. A diagonalização é feita bloco a bloco por setor de carga.

| (t = 6) | m = 0, J = 1 | m = 0,5, J = 2 |
|---|---|---|
| I(L:R) gauge/livre, qubit | 0,34 | 0,010 |
| I(L:R) gauge/livre, **inv** | **0,53** | **0,043** |
| I₃ qubit (gauge; livre) | −0,115; −0,150 | −0,096; −0,081 |
| I₃ **inv** (gauge; livre) | **−0,061; −0,001** | **−0,035; −0,003** |

| Previsão pré-registrada | Resultado |
|---|---|
| P0: entropia de região única idêntica nas duas álgebras | **Passou** (< 1e-14) |
| P1: m = 0, razão de I(L:R) maior em inv (pares do livre são férmions carregados; os do gauge, mésons neutros) | **Passou** (0,34 → 0,53) |
| P2: m = 0,5, J = 2, razão inv < 0,3 | **Passou** (0,043) |
| P3: m = 0,5, J = 2, I₃ inv do gauge negativa e abaixo da do livre | **Passou** |

**Leitura.**
1. **O efeito dinâmico do confinamento é real.** Ele aparece em quantidades de região única, que não dependem da álgebra, e sobrevive na álgebra invariante para regiões disjuntas: os pares longos caem a ~4% do livre.
2. **A álgebra invariante separa o que é carregado do que é neutro.** Os pares de quasipartículas do livre são férmions; metade da I(L:R) do livre desaparece ao proibir linhas de Wilson (1,91 → 1,00). Os mésons do gauge são neutros e perdem menos (0,64 → 0,53).
3. **Achado exploratório (não pré-registrado para m = 0): a assinatura multipartida é específica do gauge.** Na álgebra invariante, a I₃ do livre é ≈ 0, exatamente o quadro de pares de quasipartículas (emaranhamento só bipartido). A do gauge é claramente negativa, com e sem confinamento. Na álgebra de qubits, essa distinção ficava escondida (as duas eram negativas).

**Para H2.** A v0.21 mostrou que, estaticamente, o gauge só age gerando massa. Dinamicamente, ele faz algo que a massa sozinha não faz: produz **emaranhamento multipartido** (I₃ < 0) e, com confinamento, impede a formação de pares longos. Esse é o primeiro candidato a papel específico do gauge no espaço emergente que sobrevive a todas as checagens de álgebra.

## Controle com interação local (v0.23): o sinal de I₃ não é do gauge; a magnitude é

**Pergunta.** A I₃ < 0 na álgebra invariante (v0.22) é do gauge ou de qualquer interação?

**Controle.** H_V = w Σ(σ⁺σ⁻ + h.c.) + (m/2) Σ(−1)ⁿσᶻ + V Σ nᵢnᵢ₊₁ (Thirring/XXZ massivo na rede). É uma interação **local**, sem campo de gauge, com V = 1, 2, 4 e m ajustado para o **mesmo gap** de cada caso gauge. Mesmo quench exato (N = 20), mesmos blocos, álgebra invariante. Com V = 0 ele coincide exatamente com o livre.

| caso (gap) | modelo | I₃_inv (t = 6) | I₃_inv médio (t ∈ [3,6]) | I_inv(L:R) (t = 6) |
|---|---|---|---|---|
| m = 0, J = 1 (1,59) | **gauge** | **−0,061** | **−0,067** | 0,53 |
| | livre | −0,001 | −0,004 | 1,00 |
| | V = 1 / 2 / 4 | −0,015 / −0,016 / −0,002 | −0,010 / −0,009 / −0,000 | 1,56 / 1,41 / 0,63 |
| m = 0,5, J = 2 (3,45) | **gauge** | **−0,035** | **−0,038** | **0,0035** |
| | livre | −0,003 | −0,012 | 0,081 |
| | V = 1 / 2 / 4 | −0,005 / +0,002 / −0,008 | −0,010 / −0,010 / −0,006 | 0,18 / 0,12 / 0,046 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: controles com I₃_inv > −0,01 (≈ 0) | **Falhou.** V = 1 e 2 dão −0,015 e −0,016 no caso m = 0 |
| P2: \|I₃_gauge\| ≥ 3 × \|I₃_controle\| | **Passou.** 3,8–16× em t = 6; 3,7–220× na média temporal |

**Leitura.**
1. **O sinal de I₃ não é exclusivo do gauge.** A interação local também gera um pouco de emaranhamento multipartido (≈ −0,01, perto do nível do livre).
2. **A magnitude é do gauge.** A I₃ do gauge é de 4× a mais de 100× maior, com o mesmo gap e sob qualquer V testado.
3. **Exploratório: o confinamento mata pares longos muito além da interação local.** Com o mesmo gap, os pares longos do gauge são 13–53× menores que os dos controles locais. Nem a massa nem a repulsão de vizinhos reproduzem isso.

**Interpretação.** Com a lei de Gauss eliminada, o gauge em 1+1D é um potencial densidade-densidade **linear e de longo alcance**. O controle tem alcance 1. O que distingue o gauge na dinâmica é esse alcance: um potencial que cresce com a distância amarra cada partícula às demais, o que gera correlação coletiva (I₃) e impede que pares se separem.

**Para H2.** O papel específico do gauge no espaço emergente, dinamicamente, é **quantitativo e ligado ao alcance da interação**, não qualitativo. O próximo teste natural varia o alcance de um potencial linear truncado em R. Se a I₃ e a supressão de pares crescerem com R até o valor do gauge, a assinatura é "potencial confinante de longo alcance", e o gauge é o caso em que esse potencial vem de graça pela lei de Gauss.

## Potencial linear truncado, sem gauge (v0.24): cada efeito vive na escala do observável

**Identidade.** Com carga total zero, J Σₙ Lₙ² = −(J/2) Σ_{k,l} Q_k Q_l |k−l|: a energia elétrica do Schwinger é **exatamente** um potencial linear entre cargas (verificado com diferença 0,0 na matriz, N = 20).

**Família sem gauge.** H_R = hopping + massa − (J/2) Σ Q_k Q_l · min(|k−l|, R).
- R = 1 é só um termo local; R = 2 é uma interação de vizinhos; R ≥ N−1 é o Schwinger.
- Para R = 2, 3, 5, 8, 19, ajustei m para o **mesmo gap** do gauge e repeti o quench exato (N = 20) na álgebra invariante.

| R | I₃ médio t ∈ [3,6], m = 0, J = 1 | I₃ médio, m = 0,5, J = 2 | I_inv(L:R) t = 6, m = 0,5, J = 2 |
|---|---|---|---|
| 2 | −0,011 | −0,008 | 0,62 |
| 3 | −0,011 | −0,014 | 0,20 |
| 5 | −0,038 | −0,047 | 0,22 |
| 8 | −0,066 | −0,054 | 0,46 |
| 19 (= gauge) | −0,067 | −0,038 | **0,0035** |

| Previsão pré-registrada | Resultado |
|---|---|
| P0: H_{N−1} = Schwinger | **Passou** (exato) |
| P1: \|I₃\| cresce monotonicamente com R | **Falhou** no caso confinante: R = 5 e 8 *passam* do valor do gauge |
| P2: pares longos decrescem monotonicamente com R | **Falhou.** Não são monótonos, e só o alcance **completo** os suprime (≥ 20× abaixo de qualquer R < 19) |
| P3: R = 8 atinge ≥ 80% da \|I₃\| do gauge | **Passou** (97% e 143%) |

**Leitura: cada efeito é controlado pelo potencial até a escala do próprio observável.**
- **I₃** mede correlação entre três blocos que somam 6 sítios. Ela satura assim que o potencial cobre essa escala (R ≳ 5–8). Correlação multipartida local precisa de interação que alcance o tamanho da região, e não mais que isso.
- **Os pares longos** entre L e R (separados por 8 sítios) só são suprimidos quando o potencial continua crescendo além dessa separação. Um potencial truncado é como uma corda que se rompe no comprimento R: os pares escapam assim que a distância passa de R.
- A lei de Gauss em 1+1D impõe o alcance **completo** de graça. Por isso só o gauge tem os dois efeitos em todas as escalas.

**Para H2.** A assinatura dinâmica do gauge no espaço emergente é a de um **potencial confinante cujo alcance não tem corte**. Em escalas menores que R, qualquer potencial linear truncado faz o mesmo; em escalas maiores, só o gauge. Fica a pergunta para 2+1D e 3+1D, onde a lei de Gauss dá Coulomb (log ou 1/r) e não linear: o confinamento passa a depender da dinâmica não abeliana.

## Lei de escala da supressão de pares (v0.25): não é degrau, é ressonância no comprimento de ruptura

**Pergunta.** Se o potencial truncado em R age como uma corda que se rompe em R (v0.24), a supressão de pares entre blocos a distância D deveria depender só de x = R/D, como um degrau.

**Método.** Caso confinante (J = 2, gap 3,45). H_R com m ajustado ao mesmo gap, R ∈ {2, 3, 4, 5, 6, 8, 10, 12, 19}. Quench exato (N = 20). Blocos de 4 sítios com distância centro a centro D ∈ {8, 10, 12, 14}. Observável: Ī = média em t ∈ [3,6] de I_inv(A:B). Fração de escape: f = (Ī_R − Ī_gauge)/(Ī_{R=2} − Ī_gauge).

| R | Ī em D = 8 | D = 10 | D = 12 | D = 14 |
|---|---|---|---|---|
| 2 | 0,114 | 0,110 | 0,080 | 0,051 |
| 3 | 0,052 | 0,029 | 0,019 | 0,015 |
| 6 | 0,226 | 0,167 | 0,145 | 0,123 |
| 8 | **0,381** | 0,270 | 0,144 | 0,157 |
| 10 | 0,301 | **0,437** | 0,297 | 0,135 |
| 12 | 0,141 | 0,337 | **0,505** | 0,383 |
| 19 (gauge) | 0,034 | 0,005 | 0,0005 | 0,00004 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: colapso de f em R/D com forma de degrau (f ≥ 0,5 se x ≤ 0,6; f ≤ 0,2 se x ≥ 1,5) | **Falhou.** f chega a 7,45 e não é um degrau |
| P2: R_meio ∝ D | **Falhou.** f não é monótona em R, então R_meio não fica definido |

**Achado (exploratório): ressonância no comprimento de ruptura.**
- Com alcance finito, a IM entre blocos não é apenas "não suprimida": ela é **amplificada**, até 7× o caso local e 10⁴× o gauge.
- Para R ≥ 8, o **pico em D acompanha R**: 8 → 8, 10 → 10, 12 → 12.
- Interpretação: o quench cria pares que se esticam sob a tensão J até a separação R, onde a força some. Pares com energia insuficiente para escapar passam muito tempo perto do ponto de retorno (como um oscilador clássico), e os que escapam saem a partir de R. Nos dois casos o emaranhamento se concentra na distância de ruptura.
- **A escala de ruptura da corda fica impressa na geometria emergente** como um pico de conectividade.
- Há um efeito par/ímpar em R (R = 3 e 5 dão Ī baixo, com m casado ≈ 0,51 contra ≈ 0,8 para R par), ligado às células de 2 sítios da rede staggered. Ele ainda não foi isolado.

**Para H2.** O gauge (sem ruptura) suprime a conectividade em todas as distâncias. Qualquer ruptura em R cria um "atalho" de emaranhamento na escala R. Num mundo com matéria dinâmica de verdade (QCD), a corda se rompe por criação de pares. A previsão que isso sugere: **o espaço emergente de uma teoria confinante com ruptura de corda deve ter conectividade realçada na escala de ruptura** (≈ 2m/σ). É uma previsão testável no próprio Schwinger massivo, comparando cargas externas fracionárias (sem ruptura) e inteiras (com ruptura).

## Ressonância de ruptura no Schwinger real (v0.26): não existe; a da v0.25 era artefato do corte

**Pergunta.** A ressonância da v0.25 (pico de IM na distância de ruptura R) aparece no Schwinger massivo de verdade, onde a corda se rompe por criação de pares, na escala d* ≈ 2·gap/J?

**Método.** Quench exato (N = 20) sem corte artificial, nos casos (m, J) = (0,5; 2), (0,5; 1), (1; 1), (1; 0,75), com d* = 3,45, 4,97, 6,82 e 8,40. IM invariante entre duas células de 2 sítios (alinhadas à rede) a distância D = 2..14. Média em t ∈ [3,6] e também instante a instante.

| (m, J) | d* | Ī(D = 2, 4, 6, 8, 10) |
|---|---|---|
| (0,5; 2) | 3,45 | 0,249 · 0,028 · 0,009 · 0,0016 · 0,0001 |
| (0,5; 1) | 4,97 | 0,197 · 0,042 · 0,013 · 0,0031 · 0,0004 |
| (1; 1) | 6,82 | 0,188 · 0,019 · 0,003 · 0,0005 · 0,00004 |
| (1; 0,75) | 8,40 | 0,186 · 0,031 · 0,005 · 0,0009 · 0,0001 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: pico interior de Ī(D) em [4, 12] | **Falhou.** Ī(D) é monótona decrescente em todos os casos |
| P2: D_pico acompanha d* | **Falhou** (não há pico) |

Conferência: o decaimento é monótono **em todos os instantes**, então a média no tempo não escondia nenhum pico transitório.

**Leitura.**
- A ressonância da v0.25 era um **artefato do corte**: além de R a força é exatamente zero e as partículas passam a ser livres, carregando emaranhamento ao longo do caminho.
- Na ruptura física, o novo par se recombina com as extremidades em **mésons neutros e ligados**. Nenhum fragmento fica livre, e a conectividade emergente continua decaindo com a distância.
- A "nova direção" da v0.25 está **encerrada**: a escala de ruptura não fica impressa na geometria emergente do Schwinger.

**Balanço da série v0.19–v0.26 (H2 em 1+1D).**
- O que sobreviveu a todas as checagens de álgebra e de controle:
  - No vácuo, o gauge age **só gerando massa** (v0.21).
  - Na dinâmica, o potencial de Gauss de alcance ilimitado gera emaranhamento multipartido acima do que dá uma interação local (v0.22–v0.23) e suprime pares longos (v0.22, v0.24).
- O que caiu: o fator 1,5 (v0.14), a lei em m/g (v0.15), a frente rápida (v0.19) e a ressonância de ruptura (v0.25–v0.26).

## 2+1D com campo de gauge dinâmico (v0.27): a lição das cordas sobrevive fora de 1+1D

**Pergunta.** Fora de 1+1D, o campo de gauge tem dinâmica própria (plaquetes que oscilam) e a lei de Gauss deixa de ser trivial. Mesmo assim, as cordas elétricas contaminam a IM, como as linhas de Wilson implícitas da v0.19?

**Modelo.** Gauge Z₂ pura, H = −Σ Bₚ − h Σ σˣ_ℓ, com bordas lisas. Uma **dualidade exata** a leva ao Ising com campo transverso na rede dual de plaquetes (5×4): Bₚ → μˣₚ e σˣ_ℓ → μᶻₚμᶻ_q, ou seja, um μᶻ isolado é uma corda elétrica até a borda.
- Setor Z₂-par, diagonalização exata (524.288 estados), validada contra diagonalização densa (1e-15).
- Fases: desconfinada (h = 0,15), perto da transição (h ≈ h_c ≈ 0,33) e confinada (h = 0,8).

Álgebras de cada plaquete:
- **full:** os spins duais completos, o que inclui cordas μᶻμᶻ entre regiões;
- **inv:** só Bₚ, ou seja, a diagonal em μˣ.

| h | fase | IM longe/perto, full | IM longe/perto, inv | Spearman (geometria × rede), full | Spearman, **inv** |
|---|---|---|---|---|---|
| 0,15 | desconfinada | 2·10⁻⁷ | ~0 | 0,980 | **0,980** |
| 0,33 | transição | 2·10⁻³ | 1·10⁻⁷ | 0,982 | **0,983** |
| 0,8 | confinada | **0,69** | 3·10⁻⁷ | **0,574** | **0,936** |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: confinado ⇒ IM full de longo alcance (≥ 0,1) e IM inv local (≤ 1e-3) | **Passou** (0,69 contra 3·10⁻⁷) |
| P2: desconfinado ⇒ as duas decaem (≤ 1e-2) | **Passou** |
| P3: geometria inv reconstrói a rede (≥ 0,8) em todo h, e a full falha no confinado (< 0,5) | **Falhou na margem.** A inv passa (0,94–0,98); a full no confinado degrada para 0,574, mas não abaixo de 0,5 |

**Leitura.**
- No confinamento (Ising dual ordenado), o estado tem cordas elétricas condensadas. Na álgebra full elas dão IM quase constante entre plaquetes arbitrariamente distantes: é a versão 2+1D, com campo dinâmico, das linhas de Wilson implícitas da v0.19.
- A geometria emergente da álgebra full degrada justamente na fase confinada.
- Na **álgebra invariante**, a IM é local nas três fases e reconstrói a rede 2D (Spearman 0,94–0,98).
- A lição metodológica de 1+1D (v0.19–v0.20) **generaliza** para teoria de gauge com campo dinâmico: a geometria emergente de uma teoria de gauge precisa ser lida numa álgebra invariante.

**Limites.** Gauge sem matéria, rede pequena (5×4), Z₂ (não U(1)) e só o estado fundamental. O centro elétrico (fluxos de fronteira como variáveis clássicas) ainda não foi medido em 2+1D. É ali que deve aparecer o análogo de rede do termo de contato de Kabat, a ponte com G (v0.11–v0.12).

## Ponte H1 ↔ H2: centro elétrico em 2+1D (v0.28) — o emaranhamento do campo de gauge é quase todo de borda

**Pergunta.** Na v0.12, o termo de contato de Kabat (−6 por vetor) fechou S = A/4G com o G induzido. Na literatura (Donnelly-Wall 2016), ele é interpretado como **modos de borda** do campo de gauge. Qual é a peça correspondente na rede, e quanto ela pesa?

**Método.**
- Gauge Z₂ pura em 2+1D numa rede **aberta de 4×4 plaquetes, com todos os 40 elos** e lei de Gauss em todos os vértices.
- Estado fundamental exato via dual: Ising com campo de borda, 65.536 estados. Ele é mapeado em laços elétricos (bijeção μᶻ ↔ configurações σˣ).
- Para cada região A (os elos de um retângulo de plaquetes, 99 regiões), a entropia no **Hilbert estendido** se decompõe exatamente em S_ext = H_borda + S_dest:
  - H_borda é a entropia de Shannon dos fluxos elétricos nos vértices de corte, ou seja, o **centro elétrico**;
  - S_dest é a parte destilável.
- Validação: em h → 0 (toric code), S_ext = (n_corte − 1)·ln 2 **exatamente**, e toda a entropia está no centro.

| h | fase | fração de borda H_borda/S_ext (mín. / média) | H_borda vs n_corte: inclinação (R²) | intercepto de S_ext / ln 2 |
|---|---|---|---|---|
| 0,15 | desconfinada | **0,998** / 0,9998 | 0,60 (0,983) | **−0,48** |
| 0,33 | transição | 0,979 / 0,998 | 0,39 (0,953) | −0,14 |
| 0,8 | confinada | 0,995 / 0,9998 | 0,11 (0,952) | −0,06 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: no desconfinado, H_borda/S_ext ≥ 0,9 | **Passou** (0,998) |
| P2: H_borda segue lei de perímetro (R² ≥ 0,95), com inclinação decrescente em h | **Passou** |
| P3: intercepto ≈ −ln 2 no desconfinado (entropia topológica) e ≈ 0 no confinado | **Passou** (−0,48 ln 2; −0,06 ln 2) |

**Achado além do previsto.** Em **todas** as fases, 98% a 99,98% do emaranhamento do campo de gauge está no centro elétrico. A parte destilável nunca passa de 0,035 nats.

**Leitura para a ponte.**
- Na rede, o que corresponde ao termo de contato (os modos de borda) não é uma correção: para um campo de gauge sem fótons propagantes (Z₂), é **praticamente todo** o emaranhamento, com lei de perímetro (área em 2+1D).
- Isso dá suporte, no lado da rede, à leitura usada na v0.12: a contribuição de gauge à entropia de área vem da fronteira (do fluxo que atravessa a superfície), não de graus de liberdade "de bulk".
- O peso relativo borda/bulk, porém, **não é universal**: depende de haver fótons. O Z₂ não tem polarizações transversais.

**Limites.** Z₂ (discreto, sem fótons) e rede 4×4. A comparação quantitativa com −6 por vetor exige um campo com fótons (U(1)) e o limite do contínuo.

## Ponte: Z_N → U(1) em 2+1D (v0.29) — a borda continua dominando; o regime de fótons está fora do alcance exato

**Pergunta.** Na v0.28, para Z₂, 98–99,98% do emaranhamento ficou na borda (centro elétrico). Com fótons (U(1)), as polarizações transversais deveriam dar emaranhamento destilável de bulk. A fração destilável cresce quando Z_N se aproxima de U(1)?

**Método.**
- Gauge Z_N pura numa rede aberta de 3×3 plaquetes (24 elos), com H = (h/2) Σ E² − K Σ (Uₚ + Uₚ†)/2, h = 0,5 e K = 1.
- Dual exato em alturas nₚ ∈ Z_N, com E = n_esq − n_dir. Lei de Gauss sem nenhuma violação; N = 2 reproduz o Z₂ exatamente.
- N = 2..5 (até 1,95 milhão de estados). Hilbert estendido, S_ext = H_borda + S_dest, para todos os retângulos com MDR tratável.

| N | ⟨S_dest/S_ext⟩ | máx | R² (lei de perímetro de H_borda) | ⟨E²⟩ | P(\|E\| = N/2) |
|---|---|---|---|---|---|
| 2 | 0,00008 | 0,002 | 0,974 | — | — |
| 3 | 0,00011 | 0,003 | 0,973 | 0,54 | 0,62 |
| 4 | 0,00091 | 0,021 | 0,959 | 0,42 | 0,042 |
| 5 | 0,00102 | 0,019 | 0,926 | 0,37 | 0,031 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: fração destilável cresce monotonicamente com N | **Passou** no conjunto pré-registrado. Não é robusto: no conjunto comum de 25 regiões há inversão de N = 4 para 5 (0,00110 → 0,00102) |
| P2: N = 2 ≤ 0,05 e N = 5 ≥ 0,15 | **Falhou.** N = 5 dá 0,001, duas ordens abaixo |
| P3: lei de perímetro do centro (R² ≥ 0,9) | **Passou** |

**Leitura.**
- Com N = 4 ≈ 5 e P(|E| = N/2) ≈ 3%, o modelo está **quase convergido para o U(1)** nesse acoplamento. Mesmo assim, **o emaranhamento continua > 99,9% de borda**.
- Mas este é o U(1) compacto em acoplamento intermediário (g² ≈ 0,7, ⟨E²⟩ ≈ 0,37): perto do vácuo elétrico, com fótons pouco desenvolvidos numa rede 3×3.
- O regime de fótons (acoplamento fraco) exige ⟨E²⟩ grande, isto é, N muito maior, o que está fora do alcance da diagonalização exata.

**Próximo passo correto.** Ir direto ao limite de fótons com **Maxwell não compacto na rede**, que é gaussiano e exatamente solúvel em redes grandes. Para ele, Casini e Huerta mostram como calcular a entropia com centro elétrico, e o resultado pode ser comparado com o termo de contato de Kabat (−6 por vetor, v0.12) e com o modo de borda de Donnelly-Wall.

## Levantamento para T2 e checagem de novidade da I₃ (24/11/2026)

**Objetivo.** Encontrar dados experimentais publicados que permitam o primeiro teste T2 de H2, e verificar se a assinatura de I₃ (v0.22–v0.23) é nova o suficiente para ser T3.

**Experimentos com teorias de gauge em simuladores quânticos:**

| Trabalho | Plataforma | Modelo | O que foi medido | Utilidade para T2 |
|---|---|---|---|---|
| Martinez et al., *Nature* 534, 516 (2016), arXiv:1605.04570 | íons, 4 qubits, digital | **o nosso Hamiltoniano** (Gauss eliminada, ε₀ = 0, J = w, m = 0,5w), quench do vácuo nu | densidade de partículas ν(t), eco de Loschmidt, negatividade logarítmica entre metades (tomografia) | **Alta para validar a ferramenta**, baixa para H2: dados só em figuras, 4 passos de Trotter, descoerência (p ≈ 0,038/qubit/passo) e vazamento de 14–31% |
| Kokail et al., *Nature* 569, 355 (2019), arXiv:1810.03421 | íons, até 20 qubits, variacional | Schwinger (estado fundamental) | energias, gaps, variância | Estático, sem entropias de regiões |
| De et al. (2024), arXiv:2410.13815 | íons, analógico | Z₂ em 1+1D | ruptura de corda resolvida no espaço e no tempo | Modelo diferente (Z₂ com campo misto) |
| Liu et al. (2024), arXiv:2411.15443 | rede óptica | Schwinger (link quântico) | ruptura de corda, condição de ressonância | Observáveis de carga, não de entropia |
| Brydges et al., *Science* 364, 260 (2019) | íons, 20 qubits | XY de longo alcance | Rényi-2 por medições aleatórias | Mostra que S⁽²⁾ de regiões é mensurável em ~20 qubits |

**Novidade da I₃.**
- **Não é nova** a I₃ negativa após quenches: Caceffo & Alba (2023, arXiv:2305.10245) encontram I₃ < 0 em sistemas integráveis com multipletos de quasipartículas, e há também Marić & Fagotti (2022) e Berthiere (2024).
- **Não é novo** o uso de álgebras com centro para emaranhamento em gauge: Casini-Huerta-Rosabal (2014), Aoki et al. (2015), Ghosh-Soni-Trivedi (2015) e a abordagem relacional de 2025 (arXiv:2506.23459).
- **Não é nova** a medição experimental de entropias por setor de simetria: há protocolos de medições aleatórias para teorias de gauge (Bringewatt-Kunjummen-Mueller, *Quantum* 2024) e medições de emaranhamento resolvido por simetria.
- **Pode ser nova a combinação específica:** I₃ na álgebra invariante após o quench do Schwinger, comparada com uma interação local de mesmo gap e com o potencial linear truncado (v0.22–v0.24). Não encontrei esse controle na literatura. É uma novidade **estreita**: uma predição quantitativa num modelo conhecido, não uma física nova.

**Avaliação.**
1. **T2 viável agora:** reproduzir ν(t) e a negatividade de Martinez et al. (2016) com o nosso código, incluindo o modelo de erro publicado (desfase p por qubit por passo, Trotter de 4 passos). É uma **validação da ferramenta** contra dado real, não um teste de H2 em sentido forte: o próprio artigo já compara com a teoria.
2. **T2 para H2 propriamente:** não existe ainda. Nenhum experimento publicado mediu IM ou I₃ entre regiões **na álgebra invariante** de uma teoria de gauge.
3. **Caminho para T3:** formular a assinatura de I₃ como **protocolo mensurável**: Rényi-2, medições aleatórias com respeito à simetria, ~10–20 íons, erros realistas. Depois, verificar se ela sobrevive às imperfeições experimentais. Isso transforma o resultado da v0.23 numa predição que um laboratório pode testar.

## T2 de validação (v0.30): nosso código contra o experimento de Martinez et al. (2016)

**O que é.** O primeiro contato do código de H2 com dados de laboratório. Martinez et al. (*Nature* 534, 516, 2016) simularam **o nosso Hamiltoniano** de Schwinger num computador de 4 íons (J = w, m = 0,5w, ε₀ = 0), a partir do vácuo nu. Mediram a densidade de partículas ν(t) e a negatividade logarítmica E_N entre as metades.

**Protocolo.**
1. A previsão ideal foi calculada e **commitada antes de ler as figuras** (commit 3339280) e pré-registrada no Notion.
2. Os dados foram **digitalizados das figuras vetoriais (SVG)** do ar5iv. As coordenadas exatas de cada elemento vieram do próprio SVG, o ponto experimental é o centro da barra de erro, e a calibração usa as marcas dos eixos. A curva ideal do artigo começa exatamente em (0, 0), o que confirma a calibração.
3. A convenção de E_N (log₂) foi identificada pelo primeiro passo, onde a descoerência é mínima: medido 0,73, ideal 0,79 em log₂ (0,55 em ln).

| wt | ν ideal (nosso) | ν medido | Δ | E_N ideal | E_N medido |
|---|---|---|---|---|---|
| π/8 | 0,194 | 0,184 ± 0,005 | −0,011 | 0,79 | 0,73 ± 0,03 |
| π/4 | 0,449 | 0,439 ± 0,008 | −0,010 | 1,25 | 0,99 ± 0,05 |
| 3π/8 | 0,370 | 0,447 ± 0,008 | **+0,077** | 1,45 | 0,89 ± 0,04 |
| π/2 | 0,163 | 0,337 ± 0,008 | **+0,173** | 1,14 | 0,66 ± 0,04 |

| Previsão pré-registrada | Resultado |
|---|---|
| P1: nossa curva ideal = curva ideal do artigo (\|Δν\| ≤ 0,03) | **Passou** (máx. 0,0003 ao longo da curva) |
| P2: dados a \|Δν\| ≤ 0,07 da curva **ideal** | **Falhou** nos dois últimos passos (+0,077; +0,173) |
| P3: E_N medido > 0 e ≤ ideal + 0,1 | **Passou** (a descoerência só reduz o emaranhamento) |

**Exploratório (não pré-registrado): por que P2 falhou.** O experimento não implementa a evolução contínua. Ele usa **4 passos de Trotter** de π/8 e sofre **desfase**. Com a ordem Z · (pares) · (ímpares), nosso código reproduz **exatamente** os pontos de "erro de discretização" do artigo: 0,1982 / 0,4716 / 0,4235 / 0,2252. Somando o desfase publicado (p = 0,038 por qubit por passo), reproduz exatamente o "modelo de erro" deles: 0,1982 / 0,448 / 0,4333 / 0,3116. O erro quadrático médio contra os dados cai de **0,085 (ideal) para 0,015 (Trotter + desfase)**.

**Leitura.**
- O código de H2 está **validado contra dado real**: a física ideal coincide com a do artigo a 3·10⁻⁴, e o pipeline de erros reproduz o deles exatamente.
- A falha de P2 é **de protocolo, não de física**: o pré-registro comparou com a evolução ideal, mas o experimento implementa uma versão discretizada e ruidosa. Lição para os próximos T2/T3: pré-registrar **com o protocolo experimental** (passos de Trotter, ruído), não com a idealização.
- Isso **não** testa as afirmações específicas de H2 (álgebra invariante, I₃), porque nenhum experimento mediu essas grandezas. Mas dá base para o próximo passo: propor a medida de I₃ como predição T3 com protocolo realista.

## T3: a assinatura de I₃ como predição mensurável (v0.31) — não é testável na prática

**Pergunta.** A I₃ na álgebra invariante do gauge foi 4–100× maior que a de um controle com interação local de mesmo gap (v0.23). Ela sobrevive a um **protocolo experimental realista**, a ponto de virar uma predição que um laboratório possa testar? Pela regra da v0.30, o pré-registro usa o protocolo completo.

**Protocolo simulado.**
- Schwinger digital em **N = 16 qubits**, quench do vácuo nu, gauge (m = 0, J = 1) contra o controle local (V = 1, o mais negativo da v0.23), com massas ajustadas ao **mesmo gap** (1,6024).
- **Trotter Z · pares · ímpares, dt = 0,25.** É a ordem que reproduz exatamente Martinez 2016; o teste confere isso de novo.
- **Desfase** p por qubit por passo, por trajetórias (300).
- Observável: **I₃ de Rényi-2** na álgebra invariante (blocos de 2 sítios no centro).
- Medida por **sombras clássicas** de um disparo, com Haar local, estimador em lotes e desfasagem por setor no pós-processamento. O estimador foi validado sem viés: −0,2175 ± 0,046 contra −0,2197 exato.

| Previsão pré-registrada | Resultado |
|---|---|
| P1: evolução exata, \|I₃⁽²⁾_gauge\| ≥ 3 × \|I₃⁽²⁾_controle\| (média t ∈ [3,6]) | **Falhou:** −0,100 contra −0,040 (razão 2,53) |
| P2: com p = 0,01, Δ = I₃(gauge) − I₃(controle) retém ≥ 50% do valor sem ruído | **Falhou:** 42% (−0,087 → −0,036) |
| P3: com M = 10⁴ disparos, \|Δ\| ≥ 3σ | **Falhou:** σ ≈ 0,42, \|Δ\|/σ = 0,06. Para 3σ seriam **~1,2·10⁷ disparos por modelo** |

**Observações.**
- Em Rényi-2 e N = 16 o contraste é menor (2,5×) do que em von Neumann e N = 20 (até 6,5×, v0.23).
- O erro de Trotter com dt = 0,25 é grande: sem ruído, Δ(t = 4) = −0,087 contra −0,046 exato.
- O gargalo decisivo é estatístico. I₃ é uma combinação de sete entropias que se cancelam quase todas, e cada uma é estimada com variância que cresce exponencialmente com o tamanho da região.

**Leitura.** A assinatura de I₃ é **real nos modelos, mas pequena demais** para ser medida com a tecnologia atual: ~10⁷ disparos equivalem a semanas ou meses de um experimento de íons, que faz 10⁵–10⁶ medições por dia. Como T3, ela **não passa**. Isso não refuta o resultado da v0.23 (que é sobre o modelo), mas mostra que ele não vira uma predição testável.

**O que fica para H2.** Das assinaturas específicas do gauge, a que tem sinal grande é a **supressão de pares longos sob confinamento** (13–53× contra o controle local, v0.23). Ela, porém, exige regiões grandes (~6 sítios cada), o que é estatisticamente ainda mais caro. Um T3 viável teria de usar observáveis de **região única** (Rényi-2 de um bloco, barato de medir) ou correlações de **carga** (diagonais, medidas diretamente na base computacional), e não informação multipartida.

## Próximos passos

- [x] Conferir os pesos κ e o resultado N₀ + N½ − 4N₁ (verificação espectral independente; coeficiente confere com Visser 2002)
- [x] Testar a dependência de esquema (tempo próprio, Pauli-Villars, dimensional/zeta)
- [ ] Reformular o pressuposto: que setor de reguladores (compatível com "só 3 forças") fixa G > 0 e Λ ~ M_Pl? Candidatos: torres de Kaluza-Klein (Partouche & de Vaulchier 2019), espectro de cordas, QCD-like com escala de Planck
- [ ] Conferir convenções de sinal contra Adler (1982) / Frolov-Fursaev
- [x] H2: geometria a partir de teoria de gauge na rede, Z₂ com matéria (v0.4)
- [x] H2: Z₂ com matéria em redes grandes (Kitaev, exato até 2048 Majoranas) (v0.5)
- [ ] H2: Z₂ com matéria bosônica (Fradkin-Shenker) em redes grandes via PEPS/DMRG, para comparar com a matéria fermiônica
- [x] H2: U(1) com férmions: Schwinger dinâmico em 1+1D e fundo de fluxo em 2D (v0.6)
- [ ] H2: QED₃ dinâmica (U(1) compacto em 2+1D) via DMRG em cilindro
- [x] H2: SU(2) em 1+1D com férmions (v0.7)
- [x] H2: SU(3) em 1+1D (QCD 1+1D), N = 6 local e N = 8 no GitHub Actions (v0.8)
- [x] H2: teste discriminante — mesmo gap, com e sem gauge (v0.13): H0 rejeitada em U(1) e SU(2)
- [x] H2: razão assintótica das taxas em cadeias longas via DMRG (v0.14) — **≈1,5, não 2: previsão falsificada**
- [x] H2: regime confinante m ≠ 0 (v0.15) — sem assinatura de corda no vácuo; razão = f(m/g), de 1,5 a 1
- [x] H2: confinamento na dinâmica (v0.16) — S(t) separa confinamento de blindagem; a trava da IM local vem da interação gauge (também em m = 0)
- [x] H2: correlações multipartidas após o quench (v0.17) — blocos de 2 também travam; sem pico balístico no gauge; f cresce com confinamento
- [x] H2: escala do emaranhamento após o quench (v0.18) — confinamento ⇒ lei de área; blindagem ⇒ pares longos via mésons
- [x] H2: coarse-graining em escala de méson (v0.19) — sem cone balístico na IM total (picos simultâneos, artefato de corda); nas cargas, o cone é causal e igual ao livre em m = 0
- [x] H2: localidade em álgebras invariantes de gauge (v0.20) — causalidade restaurada (v ≈ 1,9); conclusão qualitativa de H2 sobrevive; números mudam por ordens de grandeza
- [x] H2: razão gauge/livre da v0.14 na álgebra invariante (v0.21) — ≈ 1,0: o fator 1,5 era artefato
- [x] H2: confinamento na dinâmica na álgebra invariante (v0.22) — efeito real; I₃ < 0 é assinatura específica do gauge (livre: I₃ ≈ 0)
- [x] H2: controle com interação local (v0.23) — sinal de I₃ não é exclusivo do gauge; magnitude 4–100× maior é; pares longos suprimidos 13–53× além da interação local
- [x] H2: potencial linear truncado em R (v0.24) — I₃ satura em R ≳ escala do bloco; supressão de pares exige alcance completo (Gauss)
- [x] H2: lei de escala da supressão de pares (v0.25) — sem colapso em degrau; ressonância: pico de IM em D ≈ R (escala de ruptura impressa na geometria)
- [x] H2: ressonância de ruptura no Schwinger real (v0.26) — não existe; a da v0.25 era artefato do corte
- [ ] H2: isolar o efeito par/ímpar em R (células staggered)
- [x] H2: 2+1D, gauge Z₂ pura (v0.27) — cordas elétricas contaminam a IM full no confinado; a álgebra inv reconstrói 2D em todas as fases
- [x] Ponte: centro elétrico em 2+1D, Z₂ (v0.28) — borda carrega 98–99,98% do emaranhamento; lei de perímetro; γ ≈ ln 2
- [x] Ponte: Z_N → U(1) (v0.29) — borda > 99,9% até N=5 (quase U(1), acoplamento intermediário); regime de fótons fora do alcance exato
- [ ] Ponte: Maxwell NÃO compacto na rede (gaussiano), 2+1D e 3+1D — entropia com centro elétrico em redes grandes; separar borda (Shannon dos fluxos) e bulk; comparar com Kabat (−6 por vetor) e Donnelly-Wall
- [ ] H2: álgebra invariante não abeliana (SU(2), SU(3): decomposição por irreps de cor) para refazer v0.7–v0.8
- [ ] H2: SU(2)/SU(3) em cadeias longas via DMRG/MPS (escala da localidade com N)
- [x] H2: localidade emergente do espectro, caso invariante por translação (v0.9); caso desordenado inconclusivo
- [x] H2: cone de luz, estrutura causal e isotropia de Lorentz emergente (v0.9)
- [x] H2: Hamiltoniano de emaranhamento geométrico + primeira lei (v0.9)
- [x] H2: passo holográfico → Einstein linearizado (rede de tensores perfeitos + Klein-Gordon em dS₂ + reconstrução do bulk) (v0.10)
- [ ] H2: Einstein linearizado em d ≥ 3 (CFT 2+1D, bulk 4D com grávitons) e redes tensoriais com laços (HaPPY completo/MERA)
- [x] Ponte H1↔H2 (v0.11): fecha para matéria; para o MP inteiro só com os modos de borda dos campos de gauge
- [x] Termo de contato de Kabat calculado pelo espectro exato no cone (holonomia) e bloco universal l = 0 na rede (v0.12)
- [ ] Maxwell completo na rede cúbica 3D com centro elétrico (Casini-Huerta) para o termo log −16/45 vs −31/45
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
- C. Cao, S. Carroll, S. Michalakis, *Space from Hilbert space*, Phys. Rev. D 95, 024031 (2017), arXiv:1606.08444
- M. Van Raamsdonk, *Building up spacetime with quantum entanglement*, Gen. Rel. Grav. 42, 2323 (2010), arXiv:1005.3035
- D. Page, W. Wootters, Phys. Rev. D 27, 2885 (1983)
- P. Calabrese, J. Cardy, J. Stat. Mech. P06002 (2004), arXiv:hep-th/0405152
- T. Jacobson, *Black hole entropy and induced gravity*, arXiv:gr-qc/9404039
- J. Cotler, G. Penington, D. Ranard, *Locality from the spectrum*, Commun. Math. Phys. 368, 1267 (2019), arXiv:1702.06142
- E. H. Lieb, D. W. Robinson, *The finite group velocity of quantum spin systems*, Commun. Math. Phys. 28, 251 (1972)
- J. J. Bisognano, E. H. Wichmann, J. Math. Phys. 16, 985 (1975); 17, 303 (1976)
- V. Eisler, I. Peschel, *Analytical results for the entanglement Hamiltonian of a free-fermion chain*, J. Phys. A 50, 284003 (2017), arXiv:1703.08126
- T. Jacobson, *Entanglement equilibrium and the Einstein equation*, PRL 116, 201101 (2016), arXiv:1505.04753
- S. Ryu, T. Takayanagi, PRL 96, 181602 (2006), arXiv:hep-th/0603001
- F. Pastawski, B. Yoshida, D. Harlow, J. Preskill, *Holographic quantum error-correcting codes* (HaPPY), JHEP 06, 149 (2015), arXiv:1503.06237
- J. de Boer, M. P. Heller, R. C. Myers, Y. Neiman, *Holographic de Sitter geometry from entanglement in CFTs*, PRL 116, 061602 (2016), arXiv:1509.00113
- J. D. Brown, M. Henneaux, Commun. Math. Phys. 104, 207 (1986)
- D. Kabat, *Black hole entropy and entropy of entanglement*, Nucl. Phys. B 453, 281 (1995), arXiv:hep-th/9503016
- W. Donnelly, A. C. Wall, *Entanglement entropy of electromagnetic edge modes*, PRL 114, 111603 (2015), arXiv:1412.1895
- L. Susskind, J. Uglum, *Black hole entropy in canonical quantum gravity and superstring theory*, PRD 50, 2700 (1994), arXiv:hep-th/9401070
- H. Casini, M. Huerta, *Entanglement entropy of a Maxwell field on the sphere*, PRD 93, 105031 (2016), arXiv:1512.06182
- W. Donnelly, A. C. Wall, *Geometric entropy and edge modes of the electromagnetic field*, PRD 94, 104053 (2016), arXiv:1506.05792
- M. Srednicki, *Entropy and area*, PRL 71, 666 (1993), arXiv:hep-th/9303048
- A. Florio, *Two-fermion negativity and confinement in the Schwinger model*, arXiv:2312.05298
- E. Fradkin, S. Shenker, *Phase diagrams of lattice gauge theories with Higgs fields*, Phys. Rev. D 19, 3682 (1979)
- A. Kitaev, *Fault-tolerant quantum computation by anyons* (código tórico), arXiv:quant-ph/9707021
- J. Schwinger, *Gauge invariance and mass II*, Phys. Rev. 128, 2425 (1962); J. Kogut, L. Susskind, Phys. Rev. D 11, 395 (1975)
- E. A. Martinez et al., *Real-time dynamics of lattice gauge theories with a few-qubit quantum computer*, Nature 534, 516 (2016), arXiv:1605.04570
- D. R. Hofstadter, Phys. Rev. B 14, 2239 (1976)
- M. C. Bañuls, K. Cichy, J. I. Cirac, K. Jansen, S. Kühn, *Efficient basis formulation for (1+1)-dimensional SU(2) lattice gauge theory*, Phys. Rev. X 7, 041046 (2017), arXiv:1707.06434
- A. Kitaev, *Anyons in an exactly solved model and beyond*, Ann. Phys. 321, 2 (2006), arXiv:cond-mat/0506438
- H. Yao, X.-L. Qi, *Entanglement entropy and entanglement spectrum of the Kitaev model*, PRL 105, 080501 (2010), arXiv:1001.1165
- S. Trebst, P. Werner, M. Troyer, K. Shtengel, C. Nayak, *Breakdown of a topological phase*, PRL 98, 070602 (2007), arXiv:cond-mat/0609048
