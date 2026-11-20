# Pré-registro L2 — certificado de otimalidade para Clifford+√T (escrito ANTES de escrever `l2.py` e de calcular qualquer estatística, 14/11/2026)

## Contexto
Weiden et al., arXiv:2609.16659 (cyclosynth), listam como aberto: uma contrapartida em Z[ζ₁₆] da relação
entre contagem de portas e expoente de denominador, que decidiria se os circuitos deles são ótimos.
Sondagem (commit d25671a): o Teorema 4.1 de Forest et al. (arXiv:1501.04944) dá, para a forma canônica,
N = 2·n_T + 3·n_√T, com N = maior expoente de denominador em β = 2cos(π/8) na matriz de Bloch.

## Dedução a testar (álgebra simples a partir da fórmula)
Custo = n_T + 3·n_√T = N − n_T ≥ N/2. Seja (t*, q*) o circuito devolvido, c* = t* + 3q*, N* = 2t* + 3q*.
Qualquer circuito mais barato (t + 3q ≤ c* − 1) com q ≥ 1 tem N ≤ 2c* − 5 = N* + 3q* − 5; os de q = 0 são
excluídos pelo ótimo em Clifford+T. Logo basta uma enumeração COMPLETA até a profundidade
D = max(0, 3q* − 5) níveis de β além de N* para provar que o circuito é ótimo.

## Dados
Resultados publicados pelos autores (github.com/mtweiden/cyclosynth, commit 0b1ffd0): `dados/u3.csv`
(SHA-256 e97df0125288e012…) e `dados/rz.csv` (d68fd20f15bfc8df…). 500 alvos × 6 precisões × 3
sintetizadores, com a sequência de portas de cada circuito. O que já foi visto antes deste registro: o
cabeçalho, as 5 primeiras linhas de cada arquivo (ε = 1e-3; q* = 2 nas duas linhas de √T) e a contagem de
linhas por grupo. Nenhuma estatística foi calculada.

## Hipóteses
- **H1 (a fórmula vale nos circuitos reais):** para toda linha, N calculado em aritmética exata a partir
  da sequência de portas é igual a 2·t + 3·q, com t e q as sílabas contadas como os autores contam
  (`scripts/_cost.py`). Em `gridsynth` e `cyclosynth_t`, q = 0 e N = 2·t.
  Falha se houver qualquer linha com N ≠ 2t + 3q. N < 2t + 3q significa circuito fora da forma canônica
  (custo superestimado) e é relatado como achado à parte.
- **H2 (o certificado é raso):** nas linhas `cyclosynth_sqrt_t` de `u3.csv` com ε ≤ 1e-6 (1500 linhas,
  nenhuma vista), a fração com D ≤ 8 (isto é, q* ≤ 4).
  - ≥ 80%: ideia VIÁVEL; segue para provar a relação entre N e a banda k do algoritmo e conferir se a
    enumeração deles é completa.
  - entre 50% e 80%: inconclusivo; relatar a distribuição de q*.
  - < 50%: ideia do certificado raso ENCERRADA.
- O limiar 8 supõe que a janela do algoritmo (w = 2 bandas de k) cobre cerca de 8 níveis de β (N ≈ 4k).
  Essa suposição NÃO está provada aqui; é o primeiro item da etapa seguinte.

## Relatar também (descritivo)
Distribuição de q* e de c*/(N*/2) por ε, nos dois arquivos; `rz.csv` com o mesmo corte de D.

## Controles (se algum falhar, nada é interpretado)
- K1: o contador de sílabas de `l2.py` reproduz as colunas t_count, q_count e cost de todas as linhas.
- K2: N exato de portas isoladas: T → 2, Q → 3, QT (= Q³) → 3, S → 0, H → 0.
- K3: mutação — inserir "HQH" no início de 200 circuitos muda N para um valor diferente de 2t + 3q
  calculado com as contagens ANTIGAS em ≥ 95% dos casos (o teste H1 tem poder).
- K4: unitariedade exata — para toda linha, M·M† = 2^h·I em Z[ζ₁₆] (a aritmética não tem erro).

## Limites declarados
- H1 confirma um teorema publicado em dados novos; não é resultado novo.
- H2 mede só a profundidade necessária. Não prova que os circuitos são ótimos: falta a completude da
  enumeração e a relação exata entre N e k.
- Uma consulta de literatura posterior a setembro de 2026 ainda não foi feita.

## Resultado (14/11/2026, rodado depois do commit 710cb48 deste pré-registro; `l2.log`, `l2_lateral.log`)
- Correção de leitura antes da rodada válida: a porta `I` (identidade, 3 linhas) não estava prevista no
  leitor; passou a ser ignorada. Nenhuma estatística tinha sido calculada.
- Controles: K1 0/18000 divergências; K2 exato; K3 200/200 mutações detectadas; K4 0/18000 falhas de
  unitariedade. Os quatro passaram.
- **H1:** N = 2t + 3q em 16371 de 18000 linhas: todas as de `cyclosynth_t`, `cyclosynth_sqrt_t` e as de
  `gridsynth` em `rz.csv`. As 1629 divergências são todas `gridsynth` em `u3.csv`, todas com N < 2t.
  Pela letra do critério H1 falha; a causa é a prevista no registro (circuito fora da forma canônica): a
  linha de base de três rotações concatenadas não é reduzida. A fórmula não foi contrariada em nenhum
  circuito em forma canônica.
- Achado lateral (descritivo, fora do pré-registro): essa linha de base conta em média 2,5 a 2,8 portas T
  acima do mínimo N/2 (1,0% a 2,9% do total; 54% dos alvos; máximo 28). Não muda a conclusão do artigo
  (9,2 contra 3,0 log2(1/ε)).
- **H2:** fração com D ≤ 8 (q* ≤ 4) em ε ≤ 1e-6: 3/1500 = 0,2% em `u3.csv` (e 3/1500 em `rz.csv`).
  Abaixo de 50% → ideia do certificado raso ENCERRADA. A mediana de q* cresce com a precisão: 3, 5, 6, 8,
  11, 13 para ε = 1e-3 … 1e-8, então D = 3q* − 5 chega a 19–34 níveis de β.
- Descritivo: mediana de c*/(N*/2) = 1,31 (ε = 1e-3) a 1,50 (ε = 1e-8) em `u3.csv`.
- Leitura: nos ótimos encontrados pelo cyclosynth o número de sílabas √T cresce com log(1/ε); uma
  varredura por expoente de denominador não certifica otimalidade em profundidade praticável. O problema
  aberto continua aberto.
