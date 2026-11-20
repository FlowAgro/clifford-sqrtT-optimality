"""Analise da L4 (ver PREREGISTRO_L4.md). Uso: python -X utf8 l4_analise.py

Reconfere cada circuito devolvido: custo pela funcao dos autores (dados/_cost_autores.py) e distancia < eps
recalculada a partir da matriz exata (l3.distancia_exata).
"""
import csv
import os
import sys
from collections import defaultdict
from statistics import mean

import l3

sys.path.insert(0, "dados")
import _cost_autores as aut

ang = {}
for r in csv.DictReader(open("dados/u3.csv", encoding="utf8")):
    if r["synthesizer"] == "cyclosynth_sqrt_t":
        ang[(r["epsilon"], r["trial"])] = (float(r["alpha"]), float(r["beta"]), float(r["gamma"]), float(r["cost"]))


def ler(arq):
    out, ruins = defaultdict(dict), 0
    for r in csv.DictReader(open(arq, encoding="utf8")):
        if not r["gates"]:
            ruins += 1
            continue
        a, b, g, pub = ang[(r["epsilon"], r["trial"])]
        d = l3.distancia_exata(r["gates"], l3.alvo(a, b, g))
        c = aut.syllable_cost(r["gates"])
        if not (d < float(r["epsilon"]) and c == float(r["cost"])):
            ruins += 1
            continue
        out[(r["janela"], r["deadline"])][int(r["trial"])] = (c, float(r["ms"]), pub)
    return out, ruins


if os.path.exists("l4_A.csv"):
    A, ruins = ler("l4_A.csv")
    print("PARTE A (eps = 1e-4) | linhas rejeitadas na reconferencia: %d" % ruins)
    a1 = A[("2", "0")]
    ig = sum(c == pub for c, ms, pub in a1.values())
    print("A1 reproducao (janela 2, sem limite): %d/%d iguais ao publicado (%.1f%%) -> %s; menores: %d; maiores: %d; ms medio %.0f"
          % (ig, len(a1), 100 * ig / len(a1), "OK" if ig >= 0.95 * len(a1) else "FALHOU",
             sum(c < pub for c, ms, pub in a1.values()), sum(c > pub for c, ms, pub in a1.values()),
             mean(ms for c, ms, pub in a1.values())))
    otimo = {int(r[1]): int(r[3]) for r in csv.reader(open("l3_B_testemunhas.csv", encoding="utf8"))}
    for jan in ("4", "6", "8"):
        x = A.get((jan, "0"), {})
        if not x:
            continue
        no_ot = sum(x[t][0] == otimo[t] for t in x)
        melhor = sum(x[t][0] < x[t][2] for t in x)
        abaixo = sum(x[t][0] < otimo[t] for t in x)
        print("A2 janela %s: n = %d | no otimo conhecido: %d (%.0f%%) | melhoraram: %d | abaixo do 'otimo' (seria erro da L3): %d | "
              "custo medio %.2f (publicado %.2f, otimo %.2f) | s por alvo: media %.1f, max %.1f"
              % (jan, len(x), no_ot, 100 * no_ot / len(x), melhor, abaixo, mean(v[0] for v in x.values()),
                 mean(v[2] for v in x.values()), mean(otimo[t] for t in x), mean(v[1] for v in x.values()) / 1000,
                 max(v[1] for v in x.values()) / 1000))
if os.path.exists("l4_B.csv"):
    B, ruins = ler("l4_B.csv")
    print("\nPARTE B (eps = 1e-7) | linhas rejeitadas na reconferencia: %d" % ruins)
    for k, x in B.items():
        pub = mean(v[2] for v in x.values())
        c = mean(v[0] for v in x.values())
        print("B janela:limite = %s:%s | n = %d | custo medio %.2f | publicado %.2f | reducao sobre o publicado %+.2f%% | "
              "melhores que o publicado: %d, piores: %d | s por alvo: media %.1f, max %.1f"
              % (k[0], k[1], len(x), c, pub, 100 * (pub - c) / pub, sum(v[0] < v[2] for v in x.values()),
                 sum(v[0] > v[2] for v in x.values()), mean(v[1] for v in x.values()) / 1000, max(v[1] for v in x.values()) / 1000))
