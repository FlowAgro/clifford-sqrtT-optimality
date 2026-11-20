"""PASSO 3 da frente MISTURA: mistura SO com T (Clifford+T) contra a mistura com sqrtT, nos MESMOS alvos.

PERGUNTA. A sonda M mostrou mistura Clifford+sqrtT a ~0,45 do deterministico publicado, mas ~0,5 e o ganho
generico de qualquer mistura (eps -> eps^2). O que seria proprio e a vantagem do sqrtT SOBREVIVER (ou crescer)
quando os dois lados misturam. No deterministico publicado, sqrtT/T = 41,41/54,57 = 0,759 em 1e-6 e
50,60/64,58 = 0,784 em 1e-7.

COMO. O mesmo enumerador da L3 (encontro no meio, sonda_A_perfil_q.preparar/candidatos) com as silabas
restritas a classe T (a = 2: HTH, SHTHs, T; eixos consecutivos distintos), isto e, a forma normal de
Matsumoto-Amano de Clifford+T vezes um Clifford. Como so ha silabas de custo 1, a tabela cresce como 3*2^n e
da para ir MUITO mais longe que no sqrtT: cobertura (a + b - 2, a mesma regra conservadora da sonda M) com
FOLGA MAIOR que a do sqrtT, o que so pode favorecer o lado T (comparacao conservadora para o sqrtT).
  Modo A: candidatos ate distancia 3e-3 (igual a sonda M), a = 18, b = 16 -> cobertura 32
          (publicado cyclosynth_t em 1e-3: media 24,84, max 28; sqrtT teve cobertura 21 contra media 18,76).
  Modo B: candidatos ate 3e-4, a = 21, b = 19 -> cobertura 38 (publicado em 1e-4: media 34,82, max 38).
  A grade de hash passa de 40r para 4r (so velocidade: a busca continua exata, celula >= 2r). Controle abaixo.
A mistura otima usa EXATAMENTE as mesmas funcoes da sonda M (pool, misturar com margem 1e-6, cota rigorosa).

CONTROLES (tem de passar, senao a rodada nao vale):
  C1 completude: em cada alvo, o menor custo de candidato a distancia < 1e-3 (modo A) / < 1e-4 (modo B) tem
     de ser <= custo publicado do cyclosynth_t nesse alvo (se der maior, falta circuito na enumeracao).
  C2 grade: no primeiro alvo, numa tabela so-T pequena (14+12, cobertura 24), a lista de candidatos com grade 4r tem de ser
     IGUAL a com grade 40r (e nao vazia).

REGRA FIXADA ANTES DE RODAR (20/11/2026), no eps efetivo decisivo (1e-6 em A, 1e-7 em B), so alvos em que as
duas misturas sao viaveis, media das razoes por alvo R = custo_mix_sqrtT / custo_mix_T:
  R <= 0,85                -> a vantagem do sqrtT SOBREVIVE a mistura (ativo possivel, segue para literatura);
  R <= 0,70                -> SINERGIA (melhor que o deterministico 0,76-0,78);
  R >  0,85                -> a vantagem some na mistura: a frente perde o valor de negocio.
  Tambem exige mistura T viavel em >= 90% dos alvos (senao: inconclusivo, cobertura insuficiente).
  Descritivo: custo_mix_T / publicado cyclosynth_t (deve dar ~0,5, controle de "mistura generica").
  O CSV sai no mesmo formato da sonda M (p com precisao total) para passar por verifica_M50.py.

AJUSTE (20/11, 6h55, teste de fumaca em 2 alvos, log FUMACA_sonda_M_A_soT_2alvos_*): o C2 em tabela 12+10
(cobertura 20) deu lista VAZIA (T-count 20 quase nunca chega a 3e-3), entao nao testava nada; passou para 14+12,
cobertura 24. Regra e cortes inalterados. A fumaca deu R ~0,76 nos 2 alvos (nao conta para o veredito).

AJUSTE 2 (20/11, 7h15, depois da 1a rodada do modo B, log sonda_M_B_soT_20261120_070520): o C2 na tabela
pequena deu lista VAZIA no raio 3e-4, e a rodada ficou INVALIDA pela regra (R = 0,751 nela, nao conta). No modo B
o C2 passa a comparar 40r contra 4r na PROPRIA tabela grande (celula 40r = 0,012 ja e fina, cabe na RAM).
Regra e cortes inalterados.

USO (WSL, venv do cyclosynth):
  python sonda_M_soT.py A --sqrt logs/sonda_M_A_20261120_042636.csv
  python sonda_M_soT.py B --sqrt logs/sonda_M_B_20261120_043031.csv
"""
import argparse
import csv
import os
import time

import numpy as np

import l3
import registro
import sonda_A_perfil_q as sa
import sonda_M_mistura as sm

AQUI = os.path.dirname(os.path.abspath(__file__))
CONF = {"A": dict(eps="1e-03", a=18, b=16, raio=3e-3, ef=["1e-05", "1e-06"], decide="1e-06"),
        "B": dict(eps="1e-04", a=21, b=19, raio=3e-4, ef=["1e-07", "1e-08"], decide="1e-07")}


def gerar_T(amax):
    """Como l3.gerar, mas so silabas de classe T (a = 2, custo 1)."""
    niv = {0: (np.array([1 + 0j]), np.array([0j]), np.array([-1], dtype=np.int8),
               np.array([-1], dtype=np.int64), np.array([-1], dtype=np.int8))}
    off, total = {0: 0}, 1
    for j in range(1, amax + 1):
        partes = []
        al, be, last, _, _ = niv[j - 1]
        for p in range(3):
            m = np.nonzero(last != p)[0]
            sa_, sb_ = l3.SIL[(p, 2)]
            partes.append((al[m] * sa_ - np.conj(be[m]) * sb_, be[m] * sa_ + np.conj(al[m]) * sb_,
                           np.full(len(m), p, dtype=np.int8), off[j - 1] + m, np.full(len(m), p * 3 + 1, dtype=np.int8)))
        niv[j] = tuple(np.concatenate([x[i] for x in partes]) for i in range(5))
        off[j] = total
        total += len(niv[j][0])
    al, be, _, par, sil = (np.concatenate([niv[j][i] for j in range(amax + 1)]) for i in range(5))
    custo = np.concatenate([np.full(len(niv[j][0]), j, dtype=np.int8) for j in range(amax + 1)])
    return al, be, par, sil, custo


def regrade(T, fator):
    """Refaz a grade de hash com celula fator*r (fator >= 2 mantem a busca exata)."""
    r = T["r"]
    hc = fator * r
    Mb = int(2 / hc) + 3
    pesos = np.array([Mb**3, Mb**2, Mb, 1], dtype=np.int64)
    chave = np.floor((T["qR"] + 1) / hc).astype(np.int64) @ pesos
    ordem = np.argsort(chave, kind="stable")
    T.update(hc=hc, pesos=pesos, KR=chave[ordem], ordem=ordem)
    return T


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("modo", choices=["A", "B"])
    ap.add_argument("--sqrt", required=True, help="CSV verificado da sonda M com sqrtT (mesmos alvos)")
    ap.add_argument("--n", type=int, default=None, help="limita o numero de alvos (teste)")
    args = ap.parse_args()
    log = registro.iniciar_log("sonda_M_%s_soT" % args.modo)
    cf = CONF[args.modo]
    cob = cf["a"] + cf["b"] - 2
    l3_gerar = l3.gerar
    l3.gerar = gerar_T
    t0 = time.time()
    T = sa.preparar(cf["a"], cf["b"], cf["raio"])
    l3.gerar = l3_gerar
    print("tabela so-T: %d palavras a esquerda, %d a direita (x24 Cliffords), cobertura %d, raio %.0e | %.0f s"
          % (len(T["al"]), T["nR"], cob, cf["raio"], time.time() - t0), flush=True)

    u3 = list(csv.DictReader(open(os.path.join(AQUI, "dados", "u3.csv"), encoding="utf8")))
    pubT = {int(r["trial"]): r for r in u3 if r["synthesizer"] == "cyclosynth_t" and r["epsilon"] == cf["eps"]}
    pubS = {int(r["trial"]): r for r in u3 if r["synthesizer"] == "cyclosynth_sqrt_t" and r["epsilon"] == cf["eps"]}
    medT = {e: np.mean([float(r["cost"]) for r in u3 if r["synthesizer"] == "cyclosynth_t" and r["epsilon"] == e])
            for e in cf["ef"]}
    sq = {}
    for r in csv.DictReader(open(args.sqrt, encoding="utf8")):
        sq[(int(r["trial"]), r["eps_ef"])] = float(r["custo_esperado"]) if r["custo_esperado"] else None
    alvos = sorted({t for t, _ in sq})[: args.n]

    # C2: grade 4r == grade 40r, numa tabela so-T PEQUENA (a = 14, b = 12; a 40r na tabela grande nao cabe na RAM)
    falhas = []
    r0 = pubS[alvos[0]]
    V0 = l3.alvo(float(r0["alpha"]), float(r0["beta"]), float(r0["gamma"]))
    if args.modo == "A":
        l3.gerar = gerar_T
        Tp, cp, rot = sa.preparar(14, 12, cf["raio"]), 24, "tabela 14+12, cobertura 24"
        l3.gerar = l3_gerar
    else:  # AJUSTE 2: em 3e-4 a tabela pequena fica vazia; a grade 40r ja e fina e cabe na tabela grande
        Tp, cp, rot = T, cob, "tabela grande, cobertura %d" % cob
    c40 = sa.candidatos(Tp, V0, cp)
    c4 = sa.candidatos(regrade(Tp, 4), V0, cp)
    ok = c40 == c4 and len(c40) > 0
    print("C2 grade (%s): 40r %d candidatos, 4r %d -> %s" % (rot, len(c40), len(c4), "IGUAIS" if ok else "DIFERENTES/VAZIA"))
    if not ok:
        falhas.append("C2")
    del Tp
    regrade(T, 4)

    csv_out = log[:-4] + ".csv"
    w = csv.writer(open(csv_out, "w", encoding="utf8", newline=""))
    w.writerow(["trial", "n_candidatos", "custo_min_candidato", "eps_ef", "custo_esperado", "cota", "suporte",
                "custo_max_suporte", "frac_sqrtT", "circuitos_e_probabilidades"])
    res = {e: [] for e in cf["ef"]}
    eps_pool = float(cf["eps"])
    for tr in alvos:
        r = pubS[tr]
        V = l3.alvo(float(r["alpha"]), float(r["beta"]), float(r["gamma"]))
        ta = time.time()
        P = sm.pool(T, V, cob, cf["raio"])
        if any(x[2] for x in P):
            falhas.append("alvo %d: candidato com sqrtT na tabela so-T" % tr)
        cmin = min((x[0] for x in P), default=None)
        cmin_eps = min((x[0] for x in P if x[4] < eps_pool ** 2), default=None)
        cpub = float(pubT[tr]["cost"])
        if cmin_eps is None or cmin_eps > cpub:
            falhas.append("C1 alvo %d: menor custo < %s = %s, publicado cyclosynth_t %g" % (tr, cf["eps"], cmin_eps, cpub))
        txt = []
        for e in cf["ef"]:
            m = sm.misturar(P, float(e)) if P else None
            if m is None:
                res[e].append((tr, None))
                txt.append("%s: inviavel" % e)
                w.writerow([tr, len(P), cmin, e, "", "", "", "", "", ""])
                continue
            c, p, cota = m
            sup = [i for i in range(len(P)) if p[i] > 1e-9]
            res[e].append((tr, c))
            txt.append("%s: %.2f (suporte %d)" % (e, c, len(sup)))
            w.writerow([tr, len(P), cmin, e, "%.4f" % c, "%.3e" % cota, len(sup), max(P[i][0] for i in sup), "0.000",
                        " ".join("%s:%r" % (P[i][5], float(p[i])) for i in sup)])
        print("alvo %d: %d candidatos so-T (min %s; min < eps %s; publicado T %g) | %s | %.0f s"
              % (tr, len(P), cmin, cmin_eps, cpub, " | ".join(txt), time.time() - ta), flush=True)

    print("\n=== RESUMO (modo %s, %d alvos, %.0f s) ===" % (args.modo, len(alvos), time.time() - t0))
    for e in cf["ef"]:
        vT = [(t, c) for t, c in res[e] if c is not None]
        par = [(c, sq.get((t, e))) for t, c in vT if sq.get((t, e))]
        Rs = [s / c for c, s in par]
        print("eps efetivo %s: mistura T viavel %d/%d | custo medio mix T %.2f (publicado T %.2f -> %.3f) | "
              "pares %d: R = mix sqrtT / mix T media %.3f mediana %.3f [min %.3f, max %.3f]"
              % (e, len(vT), len(res[e]), np.mean([c for _, c in vT]) if vT else float("nan"), medT[e],
                 (np.mean([c for _, c in vT]) / medT[e]) if vT else float("nan"), len(par),
                 np.mean(Rs) if Rs else float("nan"), np.median(Rs) if Rs else float("nan"),
                 min(Rs) if Rs else float("nan"), max(Rs) if Rs else float("nan")))
    e = cf["decide"]
    vT = [c for _, c in res[e] if c is not None]
    Rs = [sq[(t, e)] / c for t, c in res[e] if c is not None and sq.get((t, e))]
    via = len(vT) / len(res[e])
    print("\ncontroles:", "OK" if not falhas else "%d FALHAS" % len(falhas))
    for f in falhas[:20]:
        print("  ", f)
    if falhas:
        print("VEREDITO PASSO 3: rodada INVALIDA (controle falhou)")
    elif via < 0.9 or not Rs:
        print("VEREDITO PASSO 3 (eps %s): INCONCLUSIVO, mistura T viavel em %.0f%% (< 90%%): cobertura insuficiente" % (e, 100 * via))
    else:
        R = float(np.mean(Rs))
        v = "SINERGIA" if R <= 0.70 else "vantagem do sqrtT SOBREVIVE" if R <= 0.85 else "vantagem SOME na mistura"
        print("VEREDITO PASSO 3 (eps %s): R = %.3f (cortes 0,85 / 0,70), T viavel %.0f%% -> %s" % (e, R, 100 * via, v))
    print("tabela por alvo em:", csv_out)


if __name__ == "__main__":
    main()
