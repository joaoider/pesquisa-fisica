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

**N = 8 no GitHub Actions (2 704 156 estados; runner público 4 vCPU / 15 GB; 38 min).** [Execução](https://github.com/joaoider/pesquisa-fisica/actions/runs/44568396438), resultados em `results/su3_n8.json`, conferidos por `tests/h2/test_h2_resultados_pesados.py`.

| J | Energia | Singleto (S²) | Cadeia | Spearman | I(viz.)/I(d=3) |
|---|---|---|---|---|---|
| 0 | −15,849563 (= férmions livres) | 4×10⁻²² | ✔ | 0,903 | 9,7 |
| 0,5 | −14,204158 | 5×10⁻²² | ✔ | 0,966 | 24,6 |
| 2 | −11,121905 | 1×10⁻²² | ✔ | 0,966 | 169 |
| 5 | −8,687099 | 9×10⁻²³ | ✔ | 0,966 | 2132 |

Veredito: **P1, P2 e P3 confirmadas em N = 8.** Informação mútua média por distância (1…7): em J = 0 oscila (1,47; 0,098; **0,104**; 0,022; **0,028**; 0,005; **0,018**); em J = 0,5 decai (1,35; 0,13; 0,045; 0,012; 0,004; 0,001; 0,0007); em J = 5 cai a ~0 a partir da distância 3.

**Achado novo:** sem a força (J = 0), a informação mútua a distância 3 **supera** a de distância 2 (0,135 × 0,099): as duas sub-redes de férmions livres se correlacionam entre si, e a ordem por distância se inverte (Spearman 0,77). Em N = 4 o efeito chega a ligar as pontas da cadeia. **Qualquer J > 0 restaura a ordem** (Spearman ≥ 0,91). Ou seja, a matéria sozinha não define uma métrica bem-comportada, e a força de cor é o que a organiza.

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

Arquivo: `tests/test_ponte_h1_h2.py`. A entropia de emaranhamento de área (método das réplicas, heat kernel no cone) e o G induzido de H1 foram calculados **no mesmo esquema** (corte de tempo próprio).

| Predição | Resultado |
|---|---|
| **P1.** Termo da ponta do cone (1/12)(1/α − α) | **Passou.** Espectro exato (zeros de Bessel) em cones α = 1/2, 1/3, 1/4: erro ~10⁻⁵ |
| **P2.** Escalar mínimo: S/A = 1/(4G_ind) | **Passou, exatamente.** S/A = Λ²/48π = A/4G de H1 (10⁻⁸) |
| **P3.** Férmions: mesmo peso | **Passou.** c(Dirac) = 1,027 ≈ c(escalar) = 0,985 em 1+1D na rede; com 2 componentes transversas, Weyl = 1 unidade de escalar, igual ao peso em 1/G |
| **P4.** Escalar com ξ ≠ 0 | Diferença −6ξ (termo de Wald): o escalar conforme tem entropia mas **não induz G** |
| **P4.** Vetores de gauge | **Não fecha só com o bulk.** Contagem de estados +2 por vetor; 1/G dá −4. Falta −6 por vetor = termo de contato de Kabat = **modos de borda** do campo de gauge (Donnelly-Wall) |
| **P5.** Modelo Padrão inteiro | Sem modos de borda: **S_emaranhamento = 73 × A/(4G_ind)** (`xfail`). Com modos de borda: **razão = 1 exatamente** (com ou sem ν_R, ξ pequeno) |

**Leitura:** a ponte fecha, mas só com uma peça específica. As partes de matéria (escalares mínimos e férmions) têm entropia de emaranhamento exatamente igual a A/(4G_ind). Os campos das forças não: o que fecha a conta são os **modos de borda dos campos de gauge** na superfície de emaranhamento, com contribuição negativa (−6 por vetor). Para a tese, isso significa que, se a gravidade emerge das 3 forças, a entropia de Bekenstein-Hawking não é só "contagem de estados de matéria": ela inclui os graus de liberdade de borda das forças. Esse é o mesmo sinal negativo dos vetores em G (H1) e o mesmo papel da estrutura de gauge que não gera geometria sozinha (código tórico, v0.4).

**Pendente:** calcular o termo de contato −6 diretamente, com Maxwell na rede em espaço de Hilbert estendido, em vez de tomá-lo da literatura.

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
- [ ] H2: SU(2)/SU(3) em cadeias longas via DMRG/MPS (escala da localidade com N)
- [x] H2: localidade emergente do espectro, caso invariante por translação (v0.9); caso desordenado inconclusivo
- [x] H2: cone de luz, estrutura causal e isotropia de Lorentz emergente (v0.9)
- [x] H2: Hamiltoniano de emaranhamento geométrico + primeira lei (v0.9)
- [x] H2: passo holográfico → Einstein linearizado (rede de tensores perfeitos + Klein-Gordon em dS₂ + reconstrução do bulk) (v0.10)
- [ ] H2: Einstein linearizado em d ≥ 3 (CFT 2+1D, bulk 4D com grávitons) e redes tensoriais com laços (HaPPY completo/MERA)
- [x] Ponte H1↔H2 (v0.11): fecha para matéria; para o MP inteiro só com os modos de borda dos campos de gauge
- [ ] Termo de contato de Kabat calculado na rede (Maxwell com centro elétrico)
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
- E. Fradkin, S. Shenker, *Phase diagrams of lattice gauge theories with Higgs fields*, Phys. Rev. D 19, 3682 (1979)
- A. Kitaev, *Fault-tolerant quantum computation by anyons* (código tórico), arXiv:quant-ph/9707021
- J. Schwinger, *Gauge invariance and mass II*, Phys. Rev. 128, 2425 (1962); J. Kogut, L. Susskind, Phys. Rev. D 11, 395 (1975)
- E. A. Martinez et al., *Real-time dynamics of lattice gauge theories with a few-qubit quantum computer*, Nature 534, 516 (2016), arXiv:1605.04570
- D. R. Hofstadter, Phys. Rev. B 14, 2239 (1976)
- M. C. Bañuls, K. Cichy, J. I. Cirac, K. Jansen, S. Kühn, *Efficient basis formulation for (1+1)-dimensional SU(2) lattice gauge theory*, Phys. Rev. X 7, 041046 (2017), arXiv:1707.06434
- A. Kitaev, *Anyons in an exactly solved model and beyond*, Ann. Phys. 321, 2 (2006), arXiv:cond-mat/0506438
- H. Yao, X.-L. Qi, *Entanglement entropy and entanglement spectrum of the Kitaev model*, PRL 105, 080501 (2010), arXiv:1001.1165
- S. Trebst, P. Werner, M. Troyer, K. Shtengel, C. Nayak, *Breakdown of a topological phase*, PRL 98, 070602 (2007), arXiv:cond-mat/0609048
