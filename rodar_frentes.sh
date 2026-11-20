#!/usr/bin/env bash
# Roda as tres sondas em sequencia; cada uma grava seu log em l2_ciclotomico/logs/.
# Uso (no WSL):  bash ./rodar_frentes.sh [previa|completa]
#   previa   (padrao, ~40 min): teto com 5 alvos por eps, sonda A com 40 alvos em 1e-3, sonda C completa
#   completa (~3 h): sonda M em 1e-3, sonda A em 1e-4, sonda M em 1e-4, depois teto com 30 alvos por eps
#   mistura  (~30 min?): so a sonda M (1e-3 com 40 alvos, 1e-4 com 20)
set -u
MODO="${1:-previa}"
source ~/cyclosynth_l3/.venv/bin/activate
cd . || exit 1
mkdir -p logs
echo "== modo $MODO | $(date) =="

if [ "$MODO" != "completa" ]; then  # a sonda C ja rodou inteira em 19/11; na completa nao repete
  echo "== 1/3 sonda C (angulos distintos em circuitos publicos) =="
  python sonda_C_angulos.py || echo "!! sonda C falhou (ver log)"
fi

if [ "$MODO" = "mistura" ]; then
  echo "== sonda M (mistura) em 1e-3 e 1e-4 =="
  python sonda_M_mistura.py A --n 40 || echo "!! sonda M (A) falhou"
  python sonda_M_mistura.py B --n 20 || echo "!! sonda M (B) falhou"
  echo "== fim | $(date) =="; ls -1t logs | head -6; exit 0
fi

echo "== 2/3 sonda A (perfil do otimo por numero de sqrtT) =="
if [ "$MODO" = "completa" ]; then
  echo "== sonda M (mistura) em 1e-3 =="
  python sonda_M_mistura.py A --n 40 || echo "!! sonda M (A) falhou"
  python sonda_A_perfil_q.py B --n 40  || echo "!! sonda A (1e-4) falhou"
  echo "== sonda M (mistura) em 1e-4 =="
  python sonda_M_mistura.py B --n 20 || echo "!! sonda M (B) falhou"
else
  python sonda_A_perfil_q.py A --n 40  || echo "!! sonda A falhou"
fi

echo "== 3/3 teto do cyclosynth (custo x tempo em 1e-5, 1e-6, 1e-7) =="
if [ "$MODO" = "completa" ]; then
  python teto_rodar.py rodar --saida teto.csv
  python teto_rodar.py analisar --saida teto.csv
else
  python teto_rodar.py rodar --n 5 --saida teto_previa.csv
  python teto_rodar.py analisar --saida teto_previa.csv
fi
echo "== fim | $(date) | logs em ./logs =="
ls -1t logs | head -8
