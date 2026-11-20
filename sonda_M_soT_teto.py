"""SONDA M so-T COM TETO: mede quanto o corte de custo da lista (teto da enumeracao) infla o custo da mistura.

MOTIVO (auditoria do Fable, 20/11/2026, DOSSIE_FABLE_RESPOSTAS.md, achado principal). No passo 3 a mistura
sqrtT usa lista completa ate custo 21 (modo A) / 27 (modo B) e ENCOSTA nesse teto (34/40 alvos a <= 1 do teto em
1e-6; 8/20 em 1e-7). A mistura so-T usa cobertura 32 / 38 e nunca encosta. Entao o custo da mistura sqrtT pode
estar inflado, e R = 0,752 / 0,751 seriam cotas superiores. Nao da para subir a cobertura sqrtT barato; da para
fazer o contrario: cortar a lista so-T na MESMA folga do sqrtT e ver quanto o custo da mistura T sobe.

O QUE FAZ. Mesma tabela, mesmo pool e mesmo LP de sonda_M_soT.py (funcoes importadas, nada reescrito). Em cada
alvo o pool completo e calculado uma vez; depois a mistura e refeita com o pool filtrado por custo <= teto.
  Folga do sqrtT (medida nos CSVs): cobertura - deterministico publicado medio = 21 - 18,76 = 2,2 (A);
  27 - 26,24 = 0,8 (B). Teto principal so-T com a mesma folga absoluta sobre o publicado so-T (24,84 / 34,82):
    modo A: teto 27        modo B: teto 35  (folga 0,2: arredondado para BAIXO de 35,6, o corte mais duro)
  Descritivo: curva de tetos (A: 32, 30, 29, 28, 27, 26; B: 38, 37, 36, 35, 34) e a variante "por alvo":
    teto_i = det_T(alvo) + (cobertura_sqrtT - det_sqrtT(alvo)), a mesma folga absoluta alvo a alvo.

CONTROLES (tem de passar, senao a rodada nao vale):
  C0 reproducao: com o teto igual a cobertura cheia, o custo esperado de cada alvo tem de bater com o CSV da
     rodada so-T original (--soT) com diferenca <= 2e-3.
  C1 monotonia: em cada alvo o custo nao pode CAIR quando o teto desce (tolerancia 1e-6).

REGRA FIXADA ANTES DE RODAR (20/11/2026), no eps efetivo decisivo (1e-6 em A, 1e-7 em B), no teto principal,
so alvos em que a mistura T com teto e a mistura sqrtT sao viaveis:
  x      = media por alvo de custo_T(teto) / custo_T(cheio) - 1      (inflacao causada pelo teto)
  R_teto = media por alvo de custo_mix_sqrtT / custo_T(teto)         (R com os dois lados sob a mesma folga)
  x <= 1%                 -> o teto nao pesa: R = 0,75 vale como medido e "sem sinergia" esta mantido.
  x > 1% e R_teto > 0,70  -> o teto pesa, mas com os dois lados igualmente cortados ainda nao ha sinergia:
                             "sem sinergia" mantido como ESTIMATIVA (nao como medida).
  R_teto <= 0,70          -> sinergia NAO descartada: so decide subindo a cobertura sqrtT (23 no modo A).
  Se a mistura T com teto for viavel em < 90% dos alvos: INCONCLUSIVO neste teto (o corte no lado T ficou
  mais duro que o do sqrtT, que foi viavel em 100%); vale a leitura da curva.
RESSALVA escrita antes: cortar o lado T e um SUBSTITUTO de alongar o lado sqrtT. Supoe que as duas misturas
reagem ao teto de modo parecido. Nao prova nada sobre o sqrtT; mede a ordem de grandeza do vies.

USO (WSL, venv do cyclosynth; ~10 min cada modo, modo B ~1,5 GB de RAM):
  python sonda_M_soT_teto.py A --sqrt logs/sonda_M_A_20261120_042636.csv --soT logs/sonda_M_A_soT_20261120_065441.csv
  python sonda_M_soT_teto.py B --sqrt logs/sonda_M_B_20261120_043031.csv --soT logs/sonda_M_B_soT_20261120_071405.csv
Log em logs/sonda_M_<modo>_soT_teto_<data>.log e tabela em .csv do mesmo nome.
"""
import argparse
import csv
import math
import os
import time

import numpy as np

import l3
import registro
import sonda_A_perfil_q as sa
import sonda_M_mistura as sm
import sonda_M_soT as st

AQUI = os.path.dirname(os.path.abspath(__file__))
TETOS = {"A": [32, 30, 29, 28, 27, 26], "B": [38, 37, 36, 35, 34]}
PRINCIPAL = {"A": 27, "B": 35}


def ler(caminho):
    out = {}
    for r in csv.DictReader(open(caminho, encoding="utf8")):
        out[(int(r["trial"]), r["eps_ef"])] = float(r["custo_esperado"]) if r["custo_esperado"] else None
    return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("modo", choices=["A", "B"])
    ap.add_argument("--sqrt", required=True, help="CSV da sonda M com sqrtT (mesmos alvos)")
    ap.add_argument("--soT", required=True, help="CSV da rodada so-T original (controle C0)")
    ap.add_argument("--n", type=int, default=None, help="limita o numero de alvos (teste)")
    args = ap.parse_args()
    log = registro.iniciar_log("sonda_M_%s_soT_teto" % args.modo)
    cf = st.CONF[args.modo]
    cob = cf["a"] + cf["b"] - 2
    cobS = sa.CONF[args.modo]["a"] + sa.CONF[args.modo]["b"] - 2
    tetos = TETOS[args.modo]
    assert tetos[0] == cob, "o primeiro teto tem de ser a cobertura cheia"
    princ = PRINCIPAL[args.modo]

    t0 = time.time()
    l3_gerar = l3.gerar
    l3.gerar = st.gerar_T
    T = sa.preparar(cf["a"], cf["b"], cf["raio"])
    l3.gerar = l3_gerar
    st.regrade(T, 4)
    print("tabela so-T: cobertura %d, raio %.0e, cobertura sqrtT %d | tetos %s, principal %d | %.0f s"
          % (cob, cf["raio"], cobS, tetos, princ, time.time() - t0), flush=True)

    u3 = list(csv.DictReader(open(os.path.join(AQUI, "dados", "u3.csv"), encoding="utf8")))
    pubT = {int(r["trial"]): r for r in u3 if r["synthesizer"] == "cyclosynth_t" and r["epsilon"] == cf["eps"]}
    pubS = {int(r["trial"]): r for r in u3 if r["synthesizer"] == "cyclosynth_sqrt_t" and r["epsilon"] == cf["eps"]}
    sq, orig = ler(args.sqrt), ler(args.soT)
    alvos = sorted({t for t, _ in sq})[: args.n]

    csv_out = log[:-4] + ".csv"
    w = csv.writer(open(csv_out, "w", encoding="utf8", newline=""))
    w.writerow(["trial", "eps_ef", "teto", "n_candidatos", "custo_esperado", "custo_max_suporte", "custo_mix_sqrtT"])
    chaves = tetos + ["P"]
    res = {e: {k: {} for k in chaves} for e in cf["ef"]}   # res[eps][teto][trial] = custo ou None
    falhas = []
    for tr in alvos:
        r = pubS[tr]
        V = l3.alvo(float(r["alpha"]), float(r["beta"]), float(r["gamma"]))
        ta = time.time()
        P = sm.pool(T, V, cob, cf["raio"])
        teto_p = int(math.floor(float(pubT[tr]["cost"]) + cobS - float(pubS[tr]["cost"])))
        txt = []
        for e in cf["ef"]:
            ant = None
            for k in chaves:
                lim = teto_p if k == "P" else k
                Pk = [x for x in P if x[0] <= lim]
                m = sm.misturar(Pk, float(e)) if Pk else None
                c = m[0] if m else None
                res[e][k][tr] = c
                mx = max(Pk[i][0] for i in range(len(Pk)) if m[1][i] > 1e-9) if m else ""
                w.writerow([tr, e, "P%d" % lim if k == "P" else k, len(Pk), "%.4f" % c if m else "", mx,
                            sq.get((tr, e)) or ""])
                if k == cob:
                    o = orig.get((tr, e))
                    if (c is None) != (o is None) or (c is not None and abs(c - o) > 2e-3):
                        falhas.append("C0 alvo %d eps %s: agora %s, original %s" % (tr, e, c, o))
                if k != "P":
                    if ant is not None and (c is not None and c < ant - 1e-6):
                        falhas.append("C1 alvo %d eps %s teto %s: custo caiu (%.6f < %.6f)" % (tr, e, k, c, ant))
                    if ant is not None and c is None:
                        ant = float("inf")
                    elif c is not None:
                        ant = c
            cc, cp = res[e][cob][tr], res[e][princ][tr]
            txt.append("%s: cheio %s, teto %d %s" % (e, "%.2f" % cc if cc else "inviavel", princ,
                                                      "%.2f" % cp if cp else "inviavel"))
        print("alvo %d: %d candidatos (teto por alvo %d) | %s | %.0f s"
              % (tr, len(P), teto_p, " | ".join(txt), time.time() - ta), flush=True)

    print("\n=== RESUMO (modo %s, %d alvos, %.0f s) ===" % (args.modo, len(alvos), time.time() - t0))
    stat = {}
    for e in cf["ef"]:
        print("eps efetivo %s:" % e)
        for k in chaves:
            via = [t for t in alvos if res[e][k][t] is not None]
            par = [t for t in via if res[e][cob][t] is not None]
            xs = [res[e][k][t] / res[e][cob][t] - 1 for t in par]
            Rs = [sq[(t, e)] / res[e][k][t] for t in via if sq.get((t, e))]
            stat[(e, k)] = (len(via), xs, Rs)
            print("  teto %-4s viavel %2d/%d | custo medio %6.2f | inflacao x: media %5.2f%% max %5.2f%% | "
                  "R_teto (%d pares): media %.3f mediana %.3f [min %.3f, max %.3f]"
                  % ("alvo" if k == "P" else k, len(via), len(alvos),
                     np.mean([res[e][k][t] for t in via]) if via else float("nan"),
                     100 * np.mean(xs) if xs else float("nan"), 100 * max(xs) if xs else float("nan"), len(Rs),
                     np.mean(Rs) if Rs else float("nan"), np.median(Rs) if Rs else float("nan"),
                     min(Rs) if Rs else float("nan"), max(Rs) if Rs else float("nan")))

    e = cf["decide"]
    nvia, xs, Rs = stat[(e, princ)]
    print("\ncontroles:", "OK" if not falhas else "%d FALHAS" % len(falhas))
    for f in falhas[:20]:
        print("  ", f)
    if falhas:
        print("VEREDITO TETO: rodada INVALIDA (controle falhou)")
    elif nvia < 0.9 * len(alvos) or not Rs or not xs:
        print("VEREDITO TETO (eps %s, teto %d): INCONCLUSIVO, mistura T com teto viavel em %d/%d (< 90%%); ler a curva"
              % (e, princ, nvia, len(alvos)))
    else:
        x, R = float(np.mean(xs)), float(np.mean(Rs))
        v = ("sinergia NAO descartada: subir a cobertura sqrtT" if R <= 0.70 else
             "teto NAO pesa: 'sem sinergia' mantido como medido" if x <= 0.01 else
             "teto pesa, mas R_teto > 0,70: 'sem sinergia' mantido como ESTIMATIVA")
        print("VEREDITO TETO (eps %s, teto %d): x = %.2f%% (corte 1%%), R_teto = %.3f (corte 0,70), viavel %d/%d -> %s"
              % (e, princ, 100 * x, R, nvia, len(alvos), v))
    print("tabela em:", csv_out)


if __name__ == "__main__":
    main()
