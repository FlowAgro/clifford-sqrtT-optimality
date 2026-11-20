"""TETO: qual e a folga real do cyclosynth em precisao alta, e ela e alcancavel sem gastar muito mais tempo?

PERGUNTA. Um "ativo" so existe se (1) a folga de custo do modo padrao for grande o bastante e
(2) nao der para fecha-la so mexendo nos botoes que a biblioteca ja tem.

REGRAS FIXADAS ANTES DE RODAR (19/11/2026):
- Alvos NOVOS: trials 200..229 de dados/u3.csv (nunca usados em L3/L4), eps = 1e-5, 1e-6, 1e-7.
- Configuracoes (janela:limite_ms; "pad" = padrao da biblioteca; mesmos argumentos dos autores):
    pad:pad   padrao nesta maquina (so informativo; ver CORRECAO)
    pad:pad:rep  repeticao identica, para medir o ruido do limite de tempo de relogio
    4:pad     8:pad     janela maior, MESMO limite de tempo  -> "botao barato"
    2:5000    8:5000    um pouco mais de tempo
    8:60000   esforco maximo desta rodada -> serve de "melhor achado"
- Cada circuito e reconferido de forma independente: distancia exata < eps (l3.distancia_exata)
  e custo pela funcao dos autores (dados/_cost_autores.py). Linha que falha e descartada e contada.
- "Melhor achado" por alvo = menor custo entre todas as configuracoes. E um LIMITE SUPERIOR do otimo,
  entao a folga medida e um limite INFERIOR da folga verdadeira (o modo certificado dos autores
  precisaria de horizonte ~ custo em meias-unidades, inviavel em 1e-7).
- CORRECAO DE DESENHO (19/11/2026, 17h30, depois da previa de 5 alvos e ANTES da rodada completa):
  a previa mostrou que o padrao local e pior que o publicado (1e-7: 264 contra 254), porque o limite de
  tempo e de relogio e esta maquina e mais lenta que a dos autores. Logo a referencia passa a ser o
  CUSTO PUBLICADO (dados/u3.csv). Os limiares (2%, 50%, 2x) nao mudam. Alvos 200-204 entram de novo.
- Melhor achado por alvo = menor custo entre o publicado e todas as configuracoes.
- Folga(eps) = 1 - soma(melhor) / soma(publicado).
- CORTE 1 (teto): se Folga(1e-7) < 2%, nao ha ativo de custo; encerrar a frente.
- CORTE 2 (botao): se alguma configuracao com tempo medio <= 2x o tempo PUBLICADO (duration_ms) recuperar
  >= 50% da folga (recupera = (soma publicado - soma config) / (soma publicado - soma melhor)),
  a melhora ja existe na biblioteca (vira sugestao de parametro aos autores, nao invencao).
- So se passar no CORTE 1 e FALHAR no CORTE 2 ha espaco para uma heuristica nova (Etapa 1).
- Descritivo (sem corte): nos alvos em que o melhor < publicado, o melhor tem menos sqrtT e N = 2t+3q maior?
  Ruido: em quantos alvos pad:pad e pad:pad:rep diferem.

USO (no WSL, dentro de l2_ciclotomico, com o venv do cyclosynth):
  source ~/cyclosynth_l3/.venv/bin/activate
  cd .
  python teto_rodar.py rodar                 # rodada completa (~1,5 a 2,5 h; retoma de onde parou)
  python teto_rodar.py rodar --n 5           # teste rapido com 5 alvos por eps (~15 min)
  python teto_rodar.py analisar              # relatorio + veredito (pode rodar a qualquer momento)
Nao use o computador para tarefas pesadas durante a rodada: o limite de tempo e de relogio.
"""
import argparse
import csv
import os
import random
import sys
import time
from collections import defaultdict

AQUI = os.path.dirname(os.path.abspath(__file__))
EPS = ["1e-05", "1e-06", "1e-07"]
TRIALS = list(range(200, 230))
CONFS = ["pad:pad", "pad:pad:rep", "4:pad", "8:pad", "2:5000", "8:5000", "8:60000"]
CAMPOS = ["epsilon", "trial", "conf", "t_count", "q_count", "cost", "lde", "distance", "ms", "custo_publicado", "gates"]


def ler_alvos():
    alvos = {}
    for r in csv.DictReader(open(os.path.join(AQUI, "dados", "u3.csv"), encoding="utf8")):
        if r["synthesizer"] == "cyclosynth_sqrt_t":
            alvos[(r["epsilon"], int(r["trial"]))] = r
    return alvos


def rodar(args):
    import cyclosynth

    alvos = ler_alvos()
    feitos = set()
    if os.path.exists(args.saida):
        for r in csv.DictReader(open(args.saida, encoding="utf8")):
            feitos.add((r["epsilon"], int(r["trial"]), r["conf"]))
    trials = TRIALS[: args.n]
    total = len(args.eps) * len(trials) * len(CONFS)
    feito_agora, t_ini = 0, time.time()
    novo = not os.path.exists(args.saida)
    with open(args.saida, "a", encoding="utf8", newline="") as f:
        w = csv.writer(f)
        if novo:
            w.writerow(CAMPOS)
        for eps in args.eps:
            for conf in CONFS:
                jan, dl = conf.split(":")[:2]
                kw = dict(epsilon=float(eps), sqrt_t=True, near_clifford_escape=False)
                if jan != "pad":
                    kw["lde_window"] = int(jan)
                if dl != "pad":
                    kw["deadline_ms"] = int(dl)
                synth = cyclosynth.Synthesizer(**kw)
                for trial in trials:
                    feito_agora += 1
                    if (eps, trial, conf) in feitos:
                        continue
                    r = alvos.get((eps, trial))
                    if r is None:
                        print(f"  alvo ausente em u3.csv: eps={eps} trial={trial}")
                        continue
                    t0 = time.time()
                    res = synth.synthesize_u3(float(r["beta"]), float(r["alpha"]), float(r["gamma"]))
                    ms = 1000 * (time.time() - t0)
                    if res:
                        w.writerow([eps, trial, conf, res.t_count, res.q_count, res.cost, res.lde,
                                    "%.6e" % res.distance, "%.0f" % ms, r["cost"], res.gates])
                    else:
                        w.writerow([eps, trial, conf, "", "", "", "", "", "%.0f" % ms, r["cost"], ""])
                    f.flush()
                    dec = time.time() - t_ini
                    print(f"[{feito_agora}/{total}] eps={eps} conf={conf:8s} trial={trial} "
                          f"custo={res.cost if res else 'FALHOU'} ({ms/1000:.1f} s) | decorrido {dec/60:.1f} min",
                          flush=True)
    print("\nRodada terminada. Agora: python teto_rodar.py analisar")


def bootstrap_folga(pares, reps=2000, seed=1):
    """IC 95% de 1 - soma(melhor)/soma(pad), reamostrando alvos."""
    rnd = random.Random(seed)
    vals = []
    for _ in range(reps):
        amostra = [pares[rnd.randrange(len(pares))] for _ in pares]
        sp = sum(p for p, _ in amostra)
        vals.append(1 - sum(b for _, b in amostra) / sp)
    vals.sort()
    return vals[int(0.025 * reps)], vals[int(0.975 * reps)]


def analisar(args):
    import l3

    sys.path.insert(0, os.path.join(AQUI, "dados"))
    import _cost_autores as aut

    alvos = ler_alvos()
    dados = defaultdict(dict)  # eps -> {(trial, conf): (custo, ms, t, q)}
    ruins = 0
    for r in csv.DictReader(open(args.saida, encoding="utf8")):
        if not r["gates"]:
            ruins += 1
            continue
        a = alvos[(r["epsilon"], int(r["trial"]))]
        V = l3.alvo(float(a["alpha"]), float(a["beta"]), float(a["gamma"]))
        d = l3.distancia_exata(r["gates"], V)
        c = aut.syllable_cost(r["gates"])
        if not (d < float(r["epsilon"]) and c == float(r["cost"])):
            ruins += 1
            continue
        dados[r["epsilon"]][(int(r["trial"]), r["conf"])] = (c, float(r["ms"]), int(r["t_count"]), int(r["q_count"]))

    print(f"Linhas rejeitadas na reconferencia independente: {ruins}\n")
    veredito = {}
    for eps in EPS:
        D = dados.get(eps)
        if not D:
            continue
        trials = sorted({t for t, _ in D if all((t, c) in D for c in CONFS)})
        if not trials:
            print(f"eps={eps}: nenhum alvo completo ainda.\n")
            continue
        n = len(trials)
        pub = {t: float(alvos[(eps, t)]["cost"]) for t in trials}
        pub_ms = sum(float(alvos[(eps, t)]["duration_ms"]) for t in trials) / n
        pub_q = {t: int(alvos[(eps, t)]["q_count"]) for t in trials}
        pub_t = {t: int(alvos[(eps, t)]["t_count"]) for t in trials}
        melhor = {t: min([pub[t]] + [D[(t, c)][0] for c in CONFS]) for t in trials}
        s_pub, s_best = sum(pub.values()), sum(melhor.values())
        folga = 1 - s_best / s_pub
        lo, hi = bootstrap_folga([(pub[t], melhor[t]) for t in trials])
        print(f"=== eps = {eps} | n = {n} alvos completos | referencia = custo publicado ===")
        print(f"  Folga do publicado sobre o melhor achado: {100*folga:.2f}%  (IC95% bootstrap {100*lo:.2f}% a {100*hi:.2f}%)")
        print(f"  Alvos em que o publicado NAO e o melhor: {sum(pub[t] > melhor[t] for t in trials)}/{n}"
              f" | custo publicado medio {s_pub/n:.2f}, melhor medio {s_best/n:.2f}, tempo publicado {pub_ms/1000:.1f} s/alvo")
        print(f"  {'config':12s} | custo medio | vs publicado | recupera da folga | s/alvo | tempo/publicado")
        botao = []
        for c in CONFS:
            sc = sum(D[(t, c)][0] for t in trials)
            ms = sum(D[(t, c)][1] for t in trials) / n
            rec = (s_pub - sc) / (s_pub - s_best) if s_pub > s_best else float("nan")
            print(f"  {c:12s} | {sc/n:11.2f} | {100*(sc/s_pub-1):+11.2f}% | {100*rec:16.0f}% | {ms/1000:6.1f} | {ms/pub_ms:8.1f}x")
            if not c.startswith("pad:pad") and ms <= 2 * pub_ms and rec >= 0.5:
                botao.append(c)
        dif = [abs(D[(t, "pad:pad")][0] - D[(t, "pad:pad:rep")][0]) for t in trials]
        print(f"  Ruido (pad:pad x repeticao): diferem em {sum(x > 0 for x in dif)}/{n} alvos; diferenca media {sum(dif)/n:.2f}")
        menos_q = mais_N = n_mel = 0
        for t in trials:
            if melhor[t] < pub[t]:
                n_mel += 1
                cb = min(CONFS, key=lambda c: (D[(t, c)][0], D[(t, c)][1]))
                _, _, tb, qb = D[(t, cb)]
                menos_q += qb < pub_q[t]
                mais_N += (2 * tb + 3 * qb) > (2 * pub_t[t] + 3 * pub_q[t])
        if n_mel:
            print(f"  Mecanismo (alvos melhorados = {n_mel}): menos sqrtT em {menos_q}; N = 2t+3q maior em {mais_N}")
        print()
        veredito[eps] = (folga, lo, botao, n)

    print("=== VEREDITO (regras fixadas no topo do script; referencia = publicado) ===")
    if "1e-07" not in veredito:
        print("  Ainda sem dados completos em 1e-7: nao ha veredito.")
        return
    folga, lo, botao, n = veredito["1e-07"]
    if n < len(TRIALS):
        print(f"  ATENCAO: so {n} alvos em 1e-7 (rodada parcial); veredito provisorio.")
    if folga < 0.02:
        print(f"  CORTE 1 FALHOU: folga em 1e-7 = {100*folga:.2f}% < 2%. Nao ha ativo de custo. Encerrar a frente.")
    elif botao:
        print(f"  CORTE 1 passou ({100*folga:.2f}%), mas CORTE 2 mostra que {', '.join(botao)} ja recupera >= 50%")
        print("  da folga em <= 2x o tempo publicado: e ajuste de parametro da propria biblioteca, nao invencao.")
    else:
        print(f"  CORTE 1 passou ({100*folga:.2f}%; IC inferior {100*lo:.2f}%) e nenhum botao barato fecha a folga:")
        print("  ha espaco para uma heuristica nova (Etapa 1), com teto ~ essa folga. Pre-registrar antes de implementar.")


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("modo", choices=["rodar", "analisar"])
    p.add_argument("--saida", default=os.path.join(AQUI, "teto.csv"))
    p.add_argument("--n", type=int, default=len(TRIALS), help="alvos por eps (padrao 30)")
    p.add_argument("--eps", nargs="+", default=EPS, choices=EPS)
    args = p.parse_args()
    import registro
    registro.iniciar_log("teto_" + args.modo)
    rodar(args) if args.modo == "rodar" else analisar(args)


if __name__ == "__main__":
    main()
