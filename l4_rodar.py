"""L4: roda o cyclosynth (compilado no WSL, commit 0b1ffd0) variando janela de profundidade e limite de tempo.

Uso (dentro do WSL, venv de ~/cyclosynth_l3):
  python l4_rodar.py <eps> <primeiro> <ultimo> <saida.csv> <janela:deadline> [<janela:deadline> ...]
janela = "pad" usa o padrao da biblioteca; deadline = "pad" usa o padrao; 0 = sem limite de tempo.
Mesmos argumentos dos autores em scripts/gather_u3.py: sqrt_t=True, near_clifford_escape=False.
"""
import csv
import os
import sys
import time

import cyclosynth

AQUI = os.path.dirname(os.path.abspath(__file__))
eps_txt, saida = sys.argv[1], sys.argv[4]
# argv[2], argv[3] = primeiro e ultimo alvo; ou argv[2] = "lista" e argv[3] = indices separados por virgula
if sys.argv[2] == "lista":
    quais = set(int(x) for x in sys.argv[3].split(","))
else:
    quais = set(range(int(sys.argv[2]), int(sys.argv[3]) + 1))
confs = [c.split(":") for c in sys.argv[5:]]
alvos = {}
for r in csv.DictReader(open(os.path.join(AQUI, "dados", "u3.csv"), encoding="utf8")):
    if r["epsilon"] == eps_txt and r["synthesizer"] == "cyclosynth_sqrt_t" and int(r["trial"]) in quais:
        alvos[int(r["trial"])] = r
novo = not os.path.exists(saida)
with open(saida, "a", encoding="utf8", newline="") as f:
    w = csv.writer(f)
    if novo:
        w.writerow(["epsilon", "trial", "janela", "deadline", "t_count", "q_count", "cost", "distance", "ms", "custo_publicado", "gates"])
    for jan, dl in confs:
        kw = dict(epsilon=float(eps_txt), sqrt_t=True, near_clifford_escape=False)
        if jan != "pad":
            kw["lde_window"] = int(jan)
        if dl != "pad":
            kw["deadline_ms"] = int(dl)
        synth = cyclosynth.Synthesizer(**kw)
        for trial in sorted(alvos):
            r = alvos[trial]
            t0 = time.time()
            res = synth.synthesize_u3(float(r["beta"]), float(r["alpha"]), float(r["gamma"]))
            ms = 1000 * (time.time() - t0)
            w.writerow([eps_txt, trial, jan, dl, res.t_count, res.q_count, res.cost, "%.6e" % res.distance, "%.0f" % ms,
                        r["cost"], res.gates] if res else [eps_txt, trial, jan, dl, "", "", "", "", "%.0f" % ms, r["cost"], ""])
            f.flush()
