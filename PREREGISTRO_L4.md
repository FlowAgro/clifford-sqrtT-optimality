# Pré-registro L4 — mais busca no próprio cyclosynth recupera a folga? (escrito ANTES das rodadas, 14/11/2026)

## Contexto
L3 (`PREREGISTRO_L3.md`): em ε = 1e-4 há circuito mais barato que o do cyclosynth em 29 de 150 alvos, e os 51
circuitos mais baratos achados são sempre mais fundos e com menos sílabas √T. O cyclosynth expõe a janela de
profundidade (`lde_window`, padrão 2) e um limite de tempo (`deadline_ms`, padrão 600 a 2000 ms). Os autores
anotam no código que janela 3+ "dilui o limite de tempo sem recuperar custo".
Programa: cyclosynth no commit 0b1ffd0, compilado no WSL desta máquina; chamada igual à dos autores
(`sqrt_t=True, near_clifford_escape=False`), script `l4_rodar.py`.

## Já visto (sonda em alvos FORA do teste, `l4_sonda.csv`)
ε = 1e-4, alvos 490–499: padrão desta máquina dá custo médio 26,50 (publicado 26,20; 8/10 iguais); janela 2
sem limite reproduz o publicado em 10/10 (2,6 s por alvo); janelas 4, 6, 8 sem limite: 26,20 / 26,10 / 26,10
(8,9 / 18,9 / 29,4 s por alvo). ε = 1e-7, alvos 490–493 no padrão: 52,5 contra 50,5 publicado; sem limite
de tempo um alvo não terminou em vários minutos (sonda interrompida).

## Parte A — ε = 1e-4, onde o ótimo é conhecido
- A1 (controle de reprodução): janela 2 sem limite nos alvos 0–149. Esperado: custo igual ao publicado em
  ≥ 95% dos alvos. Se falhar, a Parte A não é interpretada.
- A2: nos 29 alvos subótimos da L3, janelas 4, 6 e 8 sem limite de tempo. Medida: quantos chegam ao ótimo
  conhecido (custo da testemunha da L3), quantos melhoram, tempo por alvo.
- Critério (janela 8): ≥ 80% dos 29 no ótimo → aumentar a janela resolve em 1e-4; < 30% → a janela não é a
  alavanca e seria preciso outra busca; no meio → parcial.
- Previsão: 30% a 60% dos 29 no ótimo com janela 8.

## Parte B — ε = 1e-7, onde o ótimo não é conhecido
- Alvos 0–19 de `u3.csv`. Configurações (janela:limite em ms por ramo): padrão; 2:20000; 4:20000; 4:60000.
- Medida: custo médio de cada configuração contra o custo médio PUBLICADO nos mesmos 20 alvos (a máquina
  dos autores é mais rápida; o padrão daqui serve só de referência local).
- Critério (melhor configuração): redução ≥ 3% sobre o publicado → há folga recuperável com mais busca;
  < 1% → mais tempo no algoritmo deles não recupera (ou não há folga, ou a enumeração deles não a alcança;
  este teste não distingue); no meio → inconclusivo com n = 20.
- Previsão: 1% a 3% de redução na configuração 4:60000.

## Limites declarados
- n = 20 na Parte B; o erro-padrão do custo médio é da ordem de 1%.
- Tempos medidos com 12 núcleos no WSL; não comparáveis aos do artigo.
- Reduções de custo aqui vêm do código dos autores com outros parâmetros; não são um algoritmo novo.

## Resultado (14/11/2026, depois do commit 72869d0; `l4_A.csv`, `l4_B.csv`, `l4_analise.log`)
Todos os circuitos devolvidos passaram na reconferência (custo pela função dos autores, distância < ε pela
matriz exata): 0 linhas rejeitadas.
- **A1 (controle):** janela 2 sem limite reproduz o custo publicado em 149/150 alvos (99,3%); 1 alvo saiu
  mais caro. Passou.
- **A2 (29 alvos subótimos, ε = 1e-4, sem limite de tempo):**
  | janela | no ótimo conhecido | melhoraram | custo médio | s por alvo |
  |---|---|---|---|---|
  | 4 | 8 (28%) | 9 | 27,17 | 7,9 |
  | 6 | 9 (31%) | 11 | 27,07 | 16,3 |
  | 8 | 12 (41%) | 16 | 26,86 | 26,0 |
  Referências nos mesmos 29 alvos: publicado 27,79; ótimo (L3) 26,17. Critério: 41% fica entre 30% e 80% →
  PARCIAL. Previsão (30%–60%) acertou. A janela 8 recupera 57% da folga de custo com cerca de dez vezes
  o tempo da janela 2. Nenhum circuito ficou abaixo do ótimo da L3 (validação cruzada da busca exaustiva).
- **B (20 alvos, ε = 1e-7; custo publicado médio 50,20):**
  | janela:limite (ms) | custo médio | redução sobre o publicado | s por alvo |
  |---|---|---|---|
  | padrão desta máquina | 50,75 | −1,10% | 2,4 |
  | 2:20000 | 49,20 | +1,99% | 20,9 |
  | 4:20000 | 48,90 | +2,59% | 20,9 |
  | 4:60000 | 48,70 | +2,99% | 60,9 |
  Critério: 2,99% fica abaixo do corte de 3% → INCONCLUSIVO pela letra, no limite. Previsão (1%–3%) acertou.
  A redução ainda crescia com o tempo (não saturou) e n = 20.
- Leitura: mais busca no algoritmo deles recupera parte da folga, mas devagar. A janela e o tempo não são
  uma alavanca eficiente; os circuitos mais baratos ficam fundos demais para a enumeração por níveis.

## Revisão do desenho da "busca dirigida" (14/11/2026) — ENCERRADA antes de implementar
Leitura do código-fonte do cyclosynth (commit 0b1ffd0, de 24/10/2026, posterior ao artigo):
- `src/synthesis/clifford_sqrt_t/cost_bound.rs` já contém um limite inferior certificado de custo em função
  do expoente de denominador ("cost ≥ 2N ≥ 2(2k−3)" em meias-unidades, N = expoente de Bloch reduzido), que
  é o mesmo limite custo ≥ N/2 deduzido na L2.
- `mod.rs` e `optimal.rs` já têm um modo interno de busca exaustiva com certificado de otimalidade
  (`synthesize_exhaustive_certified`, `CostCertificate`), não exposto na interface Python, além de fronteira
  ordenada por piso de custo, BKZ e poda por prefixo.
- Consequências: (1) a dedução da L2 NÃO é novidade para os autores; a afirmação anterior de que a relação
  "respondia à pergunta aberta deles" estava errada quanto à novidade — deveria ter lido o código antes;
  (2) não temos vantagem de engenharia para construir uma busca melhor que a deles; a folga medida na L3/L4
  é uma troca entre tempo e custo que eles ajustam de propósito.
- O que NÃO aparece no repositório (busca por texto): o limite de contagem 2,18·log2(1/ε) e a mistura típica
  de 22,9% de sílabas √T; a medição independente da folga (51 testemunhas por busca exaustiva própria).
