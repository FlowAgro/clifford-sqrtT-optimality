# Pré-registro L3 — o cyclosynth devolve o circuito de menor custo? (escrito ANTES de escrever `l3.py`, 14/11/2026)

## Origem (exploratório, já visto; não conta como evidência)
Contagem das formas canônicas de Clifford+√T: o número de unitárias com custo ≤ c cresce como λ^c, λ = 2,594
(raiz de x³ = 2x² + 4), o que dá o limite inferior de custo 2,18·log2(1/ε) e a mistura típica de 22,9% de
sílabas √T. Nos dados publicados do cyclosynth (`dados/u3.csv`) o ajuste é 2,40·log2(1/ε) e a fração de
sílabas √T sobe de 0,237 (ε = 1e-3) para 0,396 (ε = 1e-8). Hipótese: em precisão alta a busca deles deixa
de achar o circuito mais barato.

## Pergunta
Para os alvos publicados, existe circuito Clifford+√T dentro de ε com custo MENOR que o do cyclosynth?

## Método: busca exaustiva por encontro no meio
- Espaço: todo circuito = (palavra canônica E) · (palavra canônica D) · (Clifford), com custo(E) ≤ a e
  custo(D) ≤ b. Cobre TODOS os circuitos de custo ≤ a + b − 2 (o prefixo mais longo de custo ≤ a tem custo
  ≥ a − 2). Custo = n_T + 3·n_√T, como em `scripts/_cost.py`.
- Alvo: V = Rz(alpha)·Ry(beta)·Rz(gamma), como em `scripts/gather_u3.py`. Critério: distância diamante
  d = sqrt(1 − |tr(U†V)|²/4) < ε, estrito, sem folga.
- Todo circuito achado é reconstruído como sequência de portas e reconferido: matriz exata em Z[ζ₁₆]
  (`l2.matriz`), custo pela regra dos autores (`l2.silabas`), distância recalculada a partir da matriz exata.
  Só conta o que passa nessa reconferência.
- Não há parâmetro ajustável que mude o resultado; a busca é exaustiva.

## Dados e configuração (fixados agora)
- A: ε = 1e-3, os 500 alvos de `u3.csv`; a = 12, b = 11 (cobre custo ≤ 21; o maior c* é 22).
- B: ε = 1e-4, alvos de índice 0 a 149 de `u3.csv`; a = 16, b = 13 (cobre custo ≤ 27; certifica c* ≤ 28).
- Já visto dos dados: a distribuição de c* nessas duas precisões (para dimensionar a e b) e as médias por ε.
- Depuração permitida só nos alvos 0 a 4 de A; eles continuam contando.

## Veredito por alvo
- SUBÓTIMO: existe circuito reconferido com custo < c*.
- CERTIFICADO: nenhum circuito mais barato e c* − 1 ≤ a + b − 2.
- INDETERMINADO: nenhum mais barato achado e c* − 1 > a + b − 2.

## Critério
- Pelo menos um alvo SUBÓTIMO: achado concreto (o circuito é a prova, conferível por qualquer um).
  Relatar a fração de subótimos e a economia média.
- Zero subótimos entre os certificados: a hipótese morre nessa precisão. Continua não testada em 1e-7 e
  1e-8, onde este método não alcança.
- Previsão registrada: em A espero zero ou quase zero subótimos (a mistura bate com a contagem); em B não
  tenho previsão firme.

## Controles (se algum falhar, nada é interpretado)
- T1 convenção: a coluna `distance` dos autores é reproduzida a partir das portas deles e do meu alvo,
  com erro relativo ≤ 1e-4, em todas as linhas `cyclosynth_sqrt_t` usadas.
- T2 cobertura: para todo alvo com c* ≤ a + b − 2, a busca acha algum circuito de custo ≤ c* (o circuito
  dos autores está no espaço de busca). Uma falha indica erro na busca.
- T3 contagem: o número de palavras geradas com custo ≤ j confere com a recorrência (403660 para j = 13).
- T4 Cliffords: 24 elementos distintos a menos de fase, todos com custo 0 pela regra dos autores.

## Resultado da configuração A (14/11/2026, depois do commit 63aeefb; `l3_A.log`, `l3_A_testemunhas.csv`)
- Controles: T1 0/500 falhas; T2 0 falhas; T3 155596 palavras (confere); T4 24 Cliffords de custo 0.
- Veredito nos 500 alvos de ε = 1e-3: SUBÓTIMO 22 (4,4%), CERTIFICADO 478, INDETERMINADO 0.
- Economia nos 22: média 1,23 estado T, máximo 3. No custo total dos 500 alvos isso é 0,29%.
- Conferência independente (`l3_conferir.py`): 22/22 testemunhas confirmadas com o custo calculado pela
  função dos próprios autores (`scripts/_cost.py`, sem alteração) e a distância recalculada com 60 dígitos
  (mpmath); maior d/ε = 0,998.
- Previsão registrada ("zero ou quase zero"): 4,4% dos alvos, com efeito de 0,29% no custo.
- Leitura: em ε = 1e-3 o cyclosynth devolve o circuito de menor custo em 95,6% dos alvos; nos outros
  existe circuito mais barato dentro da mesma tolerância. Os autores não afirmam otimalidade ("we prove no
  cost-optimality"), então isto quantifica uma folga, não contradiz o artigo.

## Resultado da configuração B (14/11/2026; `l3_B.log`, `l3_B_testemunhas.csv`; 5258 s)
- Controles: T1 0/150 falhas; T2 0 falhas; T3 403660 palavras (confere); T4 24 Cliffords de custo 0.
- Veredito nos 150 alvos de ε = 1e-4: SUBÓTIMO 29 (19,3%), CERTIFICADO 120, INDETERMINADO 1 (alvo 52, c* = 29).
- Economia nos 29: 1 estado T em 17 alvos, 2 em 8, 3 em 2, 4 em 2 (média 1,62). No custo total dos 150
  alvos isso é 1,20% (47 em 3922).
- Conferência independente (`l3_conferir.py`): 29/29 testemunhas confirmadas; maior d/ε = 0,975.
- Comparação A × B: fração de subótimos 4,4% (IC 95% de Wilson 2,9%–6,6%) contra 19,3% (13,8%–26,4%);
  excesso de custo 0,29% contra 1,20%. Os intervalos não se sobrepõem: a folga cresce de 1e-3 para 1e-4.
- Não testado: ε ≤ 1e-5, onde este método não alcança. A extrapolação para 1e-8 (limite de contagem 2,18
  contra 2,40 ajustado) continua sendo hipótese.
- Exploratório, fora do registro: custo médio ótimo 18,71 (A, 500 alvos) e 25,83 (B, 150 alvos, outro
  conjunto), o que dá cerca de 2,14 por bit entre as duas precisões, perto do limite de contagem 2,18.

## Mecanismo (exploratório, fora do registro, 14/11/2026)
Nas 51 testemunhas (22 de A, 29 de B), o circuito mais barato tem SEMPRE expoente de denominador N maior
que o do cyclosynth (média +5,7 níveis de β; de +1 a +14) e SEMPRE menos sílabas √T (média −2,7 em A e
−3,0 em B), com mais sílabas T (+7,0 e +7,4). É o mecanismo previsto pela relação custo = N − n_T: os
circuitos mais baratos são mais fundos e mais ricos em T, fora da janela de profundidade da busca deles.
A interface do cyclosynth expõe `lde_window` e `deadline_ms`; não foi testado se aumentar a janela recupera
esses circuitos nem a que custo de tempo (exige compilar o programa em Rust).
