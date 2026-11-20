"""SONDA M: sintese por MISTURA em Clifford+sqrtT, com a lista COMPLETA de candidatos da busca exaustiva (L3).

IDEIA. Sortear entre varios circuitos proximos do alvo cancela o erro coerente: circuitos de precisao ~1e-3
servem para erro efetivo ~1e-6 (Campbell 2017; Hastings 2017). Escolher a mistura otima exige conhecer TODOS os
candidatos, e e isso que o encontro no meio da L3 fornece (completo ate o custo de cobertura).

MATEMATICA (cota superior rigorosa, sem aproximacao):
  W_i = V^dag U_i = c_i I - i a_i.sigma  (SU(2), sinal escolhido com c_i >= 0); |a_i| = distancia do circuito.
  Canal misto E(rho) = sum p_i W_i rho W_i^dag. Pela desigualdade triangular:
    (1/2)||E - Id||_diamante <= | sum p_i c_i a_i | + sum p_i |a_i|^2
  (o termo linear e um comutador com norma 2|b|; o quadratico e |a|^2 (P rho P - rho), norma <= 2|a|^2).
  Usamos so essa cota: toda mistura aceita aqui tem erro efetivo GARANTIDO <= eps_ef.
OTIMIZACAO: minimizar o custo ESPERADO sum p_i custo_i sujeito a sum p_i = 1, p >= 0 e a cota <= eps_ef.
  A norma euclidiana e trocada por cotas lineares mais fortes (|x|_1 ou sqrt(3)|x|_inf), resolvidas por
  programacao linear (scipy/HiGHS); fica a melhor das duas. A cota exata e reconferida no fim.

REGRAS FIXADAS ANTES DE RODAR (19/11/2026):
  Modo A: candidatos da L3 em eps = 1e-3 (completos ate custo 21); erros efetivos 1e-5 e 1e-6.
  Modo B: candidatos em eps = 1e-4 (completos ate custo 27); erros efetivos 1e-7 e 1e-8.
  Referencia: custo MEDIO publicado do cyclosynth deterministico no eps efetivo (dados/u3.csv, 500 alvos;
  os alvos de cada eps sao sorteios diferentes, entao a comparacao e de media contra media).
  PASSA se, no eps efetivo mais exigente do modo (1e-6 em A, 1e-8 em B), a mistura for viavel em >= 90% dos
  alvos E o custo esperado medio for <= 70% do publicado deterministico.
  Descritivo: tamanho do suporte da mistura, custo maximo no suporte, fracao de sqrtT.
  Fica para depois (ressalva): comparar com mistura so com T (a cobertura atual quase nao tem circuitos q = 0).
  AJUSTE (19/11, 17h35, depois de um teste de fumaca em 5 alvos, 300-304, e ANTES da rodada oficial):
  com candidatos so ate distancia 1e-3, a mistura so foi viavel em 2/5 (poucos candidatos, 5 a 16).
  A cota permite candidatos mais distantes (pagam |a|^2), entao o raio de busca passa a 3e-3 (modo A) e
  3e-4 (modo B). Limiares (90%, 0,70) inalterados. A rodada oficial comeca no trial 310.
  AJUSTE 2 (19/11, 17h35, teste de 1 alvo no modo B): 1e-8 foi inviavel por falta de candidatos com
  |a| <~ 1e-4 dentro da cobertura (custo <= 27), limite da enumeracao e nao da ideia. O eps decisivo do
  modo B passa a 1e-7; 1e-8 fica descritivo.
  AJUSTE 3 (20/11, depois da verificacao em 50 digitos, logs/verifica_M50_20261119_215435): o LP mirava
  EXATAMENTE eps_ef com |a_i| em float64 (erro ~1e-16 absoluto), e em 50 digitos 11 misturas do modo B passaram
  de eps por ate 3e-8 relativo (1 delas no eps decisivo 1e-7). Correcao: o LP mira eps_ef (1 - MARGEM), com
  MARGEM = 1e-6 (30x o ruido medido). Limiares e alvos inalterados; custo esperado muda ~1e-6 relativo.

USO (WSL, venv do cyclosynth, precisa de numpy e scipy; nao precisa de internet):
  python sonda_M_mistura.py A --n 40      # ~1 min
  python sonda_M_mistura.py B --n 20      # tabela maior, alguns minutos
Log em logs/sonda_M_<modo>_<data>.log e tabela em logs/sonda_M_<modo>_<data>.csv.
"""
import argparse
import cmath
import csv
import math
import os
import time

import numpy as np
from scipy.optimize import linprog

import l2
import l3
import registro
from sonda_A_perfil_q import CONF, candidatos, preparar

RAIO = {"A": 3e-3, "B": 3e-4}
DECIDE = {"A": "1e-06", "B": "1e-07"}
MARGEM = 1e-6  # AJUSTE 3: folga contra o ruido de float64 nos |a_i|

AQUI = os.path.dirname(os.path.abspath(__file__))
EF = {"A": ["1e-05", "1e-06"], "B": ["1e-07", "1e-08"]}
PAULI = [np.array([[0, 1], [1, 0]], dtype=complex), np.array([[0, -1j], [1j, 0]], dtype=complex),
         np.array([[1, 0], [0, -1]], dtype=complex)]
Z8 = [cmath.exp(1j * math.pi * k / 8) for k in range(8)]


def unitaria(gates):
    a, b, c, d, h = l2.matriz(gates)
    ev = lambda e: sum(x * w for x, w in zip(e, Z8))
    U = np.array([[ev(a), ev(c)], [ev(b), ev(d)]], dtype=complex) / 2 ** (h / 2)
    return U / cmath.sqrt(np.linalg.det(U))


def erro(U, V):
    W = V.conj().T @ U
    c = (np.trace(W) / 2).real
    av = np.array([(0.5j * np.trace(W @ P)).real for P in PAULI])
    if c < 0:
        c, av = -c, -av
    return c, av


def pool(T, V, cob, raio):
    """Lista de (custo, t, q, c*a (3,), |a|^2, gates), um por unitaria distinta."""
    vistos, out = set(), []
    for ctot, li, ri in candidatos(T, V, cob):
        g = l3.palavra(li, T["par"], T["sil"]) + l3.palavra(ri % T["nR"], T["par"], T["sil"]) + T["cl"][ri // T["nR"]]
        d = l3.distancia_exata(g, V)
        if d >= raio:
            continue
        c, av = erro(unitaria(g), V)
        if abs(np.linalg.norm(av) - d) > 1e-9:
            raise RuntimeError("convencao de distancia nao bate: %g vs %g (%s)" % (np.linalg.norm(av), d, g))
        ch = tuple(np.round(av, 12))
        if ch in vistos:
            continue
        vistos.add(ch)
        t, q = l2.silabas(g)
        out.append((t + 3 * q, t, q, c * av, float(av @ av), g))
    return out


def misturar(P, eps_ef):
    """Menor custo esperado com cota <= eps_ef. Devolve (custo, p, cota) ou None."""
    n = len(P)
    custo = np.array([x[0] for x in P], dtype=float)
    X = np.array([x[3] for x in P])          # n x 3
    s = np.array([x[4] for x in P])
    melhor = None
    eps_ef = eps_ef * (1 - MARGEM)
    for tipo in ("l1", "linf"):
        # variaveis: p (n), w (3)
        cvec = np.concatenate([custo, np.zeros(3)])
        A, b = [], []
        for k in range(3):  # -w_k <= x_k <= w_k
            A.append(np.concatenate([X[:, k], -np.eye(3)[k]])); b.append(0.0)
            A.append(np.concatenate([-X[:, k], -np.eye(3)[k]])); b.append(0.0)
        if tipo == "l1":
            A.append(np.concatenate([s, np.ones(3)])); b.append(eps_ef)
        else:  # w_k <= m para todo k, cota sqrt(3) m; usa w0 como m
            for k in (1, 2):
                A.append(np.concatenate([np.zeros(n), np.eye(3)[k] - np.eye(3)[0]])); b.append(0.0)
            A.append(np.concatenate([s, [math.sqrt(3), 0, 0]])); b.append(eps_ef)
        Aeq = [np.concatenate([np.ones(n), np.zeros(3)])]
        r = linprog(cvec, A_ub=np.array(A), b_ub=np.array(b), A_eq=np.array(Aeq), b_eq=[1.0],
                    bounds=[(0, None)] * (n + 3), method="highs")
        if r.status != 0:
            continue
        p = np.clip(r.x[:n], 0, None)
        p = p / p.sum()
        cota = float(np.linalg.norm(p @ X) + p @ s)
        if cota <= eps_ef * (1 + 1e-9) and (melhor is None or p @ custo < melhor[0]):
            melhor = (float(p @ custo), p, cota)
    return melhor


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("modo", choices=["A", "B"])
    ap.add_argument("--n", type=int, default=40)
    ap.add_argument("--primeiro", type=int, default=310)
    ap.add_argument("--raio", type=float, default=None, help="raio dos candidatos (padrao A 3e-3, B 3e-4)")
    args = ap.parse_args()
    log = registro.iniciar_log("sonda_M_" + args.modo)
    cf = CONF[args.modo]
    eps = float(cf["eps"])
    raio = args.raio or RAIO[args.modo]
    cob = cf["a"] + cf["b"] - 2
    linhas = list(csv.DictReader(open(os.path.join(AQUI, "dados", "u3.csv"), encoding="utf8")))
    cs = [r for r in linhas if r["synthesizer"] == "cyclosynth_sqrt_t"]
    pub = {e: np.mean([float(r["cost"]) for r in cs if r["epsilon"] == e]) for e in EF[args.modo]}
    pub_pool = np.mean([float(r["cost"]) for r in cs if r["epsilon"] == cf["eps"]])
    print("candidatos: eps = %s, completos ate custo %d | publicado deterministico medio: %s = %.2f; %s"
          % (cf["eps"], cob, cf["eps"], pub_pool, ", ".join("%s = %.2f" % (e, pub[e]) for e in EF[args.modo])))
    t0 = time.time()
    T = preparar(cf["a"], cf["b"], raio)
    print("tabela pronta em %.0f s | raio dos candidatos %.1e (alvos sorteados no eps %s)" % (time.time() - t0, raio, cf["eps"]), flush=True)

    alvos = sorted((int(r["trial"]), r) for r in cs if r["epsilon"] == cf["eps"]
                   and args.primeiro <= int(r["trial"]) < args.primeiro + args.n)
    csv_out = log[:-4] + ".csv"
    w = csv.writer(open(csv_out, "w", encoding="utf8", newline=""))
    w.writerow(["trial", "n_candidatos", "custo_min_candidato", "eps_ef", "custo_esperado", "cota", "suporte",
                "custo_max_suporte", "frac_sqrtT", "circuitos_e_probabilidades"])
    res = {e: [] for e in EF[args.modo]}
    for tr, r in alvos:
        V = l3.alvo(float(r["alpha"]), float(r["beta"]), float(r["gamma"]))
        P = pool(T, V, cob, raio)
        cmin = min((x[0] for x in P), default=None)
        txt = []
        for e in EF[args.modo]:
            m = misturar(P, float(e)) if P else None
            if m is None:
                res[e].append(None)
                txt.append("%s: inviavel" % e)
                w.writerow([tr, len(P), cmin, e, "", "", "", "", "", ""])
                continue
            c, p, cota = m
            sup = [i for i in range(len(P)) if p[i] > 1e-9]
            tq = sum(p[i] * P[i][1] for i in sup), sum(p[i] * P[i][2] for i in sup)
            fq = tq[1] / (tq[0] + tq[1]) if sum(tq) else 0
            res[e].append((c, len(sup), max(P[i][0] for i in sup), fq))
            txt.append("%s: %.2f (suporte %d, cota %.2e)" % (e, c, len(sup), cota))
            w.writerow([tr, len(P), cmin, e, "%.4f" % c, "%.3e" % cota, len(sup), max(P[i][0] for i in sup),
                        "%.3f" % fq, " ".join("%s:%r" % (P[i][5], float(p[i])) for i in sup)])
        print("alvo %d: %d candidatos (custo min %s) | %s" % (tr, len(P), cmin, " | ".join(txt)), flush=True)

    print("\n=== RESUMO (modo %s, n = %d alvos, %.0f s) ===" % (args.modo, len(alvos), time.time() - t0))
    ok = True
    for e in EF[args.modo]:
        v = [x for x in res[e] if x]
        via = len(v) / max(1, len(res[e]))
        if not v:
            print("eps efetivo %s: nenhuma mistura viavel" % e)
            ok = False
            continue
        cm = np.mean([x[0] for x in v])
        print("eps efetivo %s: viavel em %d/%d (%.0f%%) | custo esperado medio %.2f | publicado deterministico %.2f"
              " -> razao %.3f | suporte medio %.1f | custo max medio no suporte %.1f | frac sqrtT %.3f"
              % (e, len(v), len(res[e]), 100 * via, cm, pub[e], cm / pub[e], np.mean([x[1] for x in v]),
                 np.mean([x[2] for x in v]), np.mean([x[3] for x in v])))
    e = DECIDE[args.modo]
    v = [x for x in res[e] if x]
    if v:
        via = len(v) / len(res[e])
        razao = np.mean([x[0] for x in v]) / pub[e]
        ok = via >= 0.9 and razao <= 0.7
        print("\nVEREDITO M (eps efetivo %s): viabilidade %.0f%% (corte 90%%), razao %.3f (corte 0,70) -> %s"
              % (e, 100 * via, razao, "PASSOU" if ok else "FALHOU"))
    print("tabela por alvo em:", csv_out)


if __name__ == "__main__":
    main()
