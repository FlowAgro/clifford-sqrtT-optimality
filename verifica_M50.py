"""VERIFICACAO EM 50 DIGITOS das cotas da sonda M (complemento de verifica_M.py).

Por que: em float64, a cota inferior da distancia diamante saiu ate 3e-8 (relativo) ACIMA de eps e as vezes acima
da cota superior, o que e impossivel na matematica: e ruido de ~1e-16 absoluto. Aqui as duas cotas sao
calculadas com mpmath a 50 digitos, sem tolerancia de ruido.

Como (tudo que decide e feito em 50 digitos):
  - U_i: entradas EXATAS em Z[zeta16] (verifica_M.exata, inteiros), avaliadas com zeta16 = exp(i pi/8) a 50 digitos.
  - V = Rz(alpha) Ry(beta) Rz(gamma), com alpha, beta, gamma lidos de dados/u3.csv (o alvo e esse numero decimal).
  - p_i lidos dos CSVs (precisao total, %r) e renormalizados a 50 digitos.
  - Cota SUPERIOR: |sum p_i c_i a_i| + sum p_i |a_i|^2.
  - Cota INFERIOR: (1/2)||(E - Id) x I (psi psi^dag)||_1 para um psi de qubit + ancila. O psi e ESCOLHIDO por busca em
    float (Nelder-Mead, estado de Bell + 8 partidas); qualquer psi normalizado da cota inferior valida, entao a
    escolha nao precisa de precisao, so a avaliacao (feita a 50 digitos, autovalores com mpmath.eighe).

REGRA FIXADA ANTES DE RODAR (19/11/2026, ~20h):
  (a) Em TODAS as misturas viaveis dos CSVs sonda_M_A_20261119_192920 e sonda_M_B_20261119_193313:
      cota inferior <= cota superior (diferenca acima de 1e-40 absoluto = a cota esta ERRADA = falha).
  (b) Cota superior <= eps_ef (1 + 1e-9), a tolerancia do proprio programa linear, fixada em sonda_M_mistura.py
      antes das rodadas. Conta e imprime, a parte, quantas passam de eps estritamente e por quanto.
  VERIFICADA se (a) e (b) valem em todas. Descritivo: aperto = inferior/superior.

REEXECUCAO (20/11): depois do AJUSTE 3 da sonda M (margem 1e-6 no LP), a MESMA regra (a)+(b) e aplicada aos CSVs
novos, passados na linha de comando. A 1a rodada (verifica_M50_20261119_215435) deu NAO VERIFICADA no modo B.

USO (Python do Windows, que tem mpmath; ou qualquer Python com numpy, scipy e mpmath):
  python verifica_M50.py [csv_A csv_B]   # nomes em logs/; sem argumento, os CSVs de 19/11
"""
import csv
import os
import sys

import mpmath as mp
import numpy as np
from scipy.optimize import minimize

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)
import registro  # noqa: E402
import verifica_M as v  # noqa: E402  (so a aritmetica exata em Z[zeta16] e as funcoes em float para escolher psi)

mp.mp.dps = 50
CSVS = [("sonda_M_A_20261119_192920.csv", "1e-03"), ("sonda_M_B_20261119_193313.csv", "1e-04")]
I2 = mp.matrix([[1, 0], [0, 1]])
SIG = [mp.matrix([[0, 1], [1, 0]]), mp.matrix([[0, -1j], [1j, 0]]), mp.matrix([[1, 0], [0, -1]])]
Z16 = [mp.expjpi(mp.mpf(k) / 8) for k in range(8)]


def dag(A):
    return A.transpose_conj()


def tr(A):
    return sum(A[i, i] for i in range(A.rows))


def kron2(A):  # A x I_2
    K = mp.zeros(4, 4)
    for i in range(2):
        for j in range(2):
            for k in range(2):
                K[2 * i + k, 2 * j + k] = A[i, j]
    return K


def u50(gates):
    M, h = v.exata(gates)
    s = mp.sqrt(2) ** h
    return mp.matrix([[mp.fsum(M[i][j][k] * Z16[k] for k in range(8)) / s for j in range(2)] for i in range(2)])


def alvo50(al, be, ga):
    rz = lambda t: mp.matrix([[mp.exp(-0.5j * t), 0], [0, mp.exp(0.5j * t)]])
    ry = mp.matrix([[mp.cos(be / 2), -mp.sin(be / 2)], [mp.sin(be / 2), mp.cos(be / 2)]])
    return rz(al) * ry * rz(ga)


def ca50(U, V):
    W = dag(V) * U
    W = W / mp.sqrt(mp.det(W))
    c = mp.re(tr(W) / 2)
    a = [mp.re(0.5j * tr(W * P)) for P in SIG]
    return (c, a) if c >= 0 else (-c, [-x for x in a])


def kraus50(c, a):
    return c * I2 - 1j * (a[0] * SIG[0] + a[1] * SIG[1] + a[2] * SIG[2])


def f50(Ws, p, x):
    psi = mp.matrix([mp.mpc(x[k], x[k + 4]) for k in range(4)])
    psi = psi / mp.norm(psi)
    rho = psi * dag(psi)
    D = -rho
    for W, pi in zip(Ws, p):
        K = kron2(W)
        D = D + pi * (K * rho * dag(K))
    D = (D + dag(D)) / 2
    ev = mp.eighe(D, eigvals_only=True)
    return mp.fsum(abs(e) for e in ev) / 2


def melhor_psi(Wf, pf, rng):
    xs = [np.array([1, 0, 0, 1, 0, 0, 0, 0], float)] + [rng.normal(size=8) for _ in range(8)]
    melhor, xm = -1.0, xs[0]
    for x0 in xs:
        r = minimize(lambda x: -v.f_diamante(Wf, pf, x), x0, method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-10, fatol=1e-18))
        for val, x in ((-r.fun, r.x), (v.f_diamante(Wf, pf, x0), x0)):
            if val > melhor:
                melhor, xm = val, x
    return xm


def main():
    registro.iniciar_log("verifica_M50")
    rng = np.random.default_rng(20261119)
    u3 = {(r["epsilon"], int(r["trial"])): r for r in csv.DictReader(open(os.path.join(AQUI, "dados", "u3.csv"), encoding="utf8"))
          if r["synthesizer"] == "cyclosynth_sqrt_t"}
    falhas, aperto, acima = [], {}, {}
    lista = [(n, "1e-03" if "_M_A_" in n else "1e-04") for n in sys.argv[1:]] or CSVS
    for nome, pe in lista:
        print("\n== %s ==" % nome)
        for r in csv.DictReader(open(os.path.join(AQUI, "logs", nome), encoding="utf8")):
            if not r["custo_esperado"]:
                continue
            tr_, es = int(r["trial"]), r["eps_ef"]
            e = mp.mpf(es)
            ref = u3[(pe, tr_)]
            V = alvo50(mp.mpf(ref["alpha"]), mp.mpf(ref["beta"]), mp.mpf(ref["gamma"]))
            Vf = v.alvo(float(ref["alpha"]), float(ref["beta"]), float(ref["gamma"]))
            itens = [s.rsplit(":", 1) for s in r["circuitos_e_probabilidades"].split()]
            p = [mp.mpf(x) for _, x in itens]
            sp = mp.fsum(p)
            p = [x / sp for x in p]
            Ws, Wf, lin, quad = [], [], [mp.mpf(0)] * 3, mp.mpf(0)
            for (g, _), pi in zip(itens, p):
                c, a = ca50(u50(g), V)
                Ws.append(kraus50(c, a))
                lin = [l + pi * c * x for l, x in zip(lin, a)]
                quad += pi * mp.fsum(x * x for x in a)
                cf, af = v.ca(v.numerica(*v.exata(g)), Vf)
                Wf.append(v.kraus(cf, af))
            sup = mp.sqrt(mp.fsum(x * x for x in lin)) + quad
            inf = f50(Ws, p, melhor_psi(Wf, np.array([float(x) for x in p]), rng))
            msg = []
            if inf - sup > mp.mpf("1e-40"):
                msg.append("INFERIOR > SUPERIOR por %s" % mp.nstr(inf - sup, 5))
            if sup > e * (1 + mp.mpf("1e-9")):
                msg.append("superior passa de eps(1+1e-9)")
            falhas += ["%s alvo %d eps %s: %s" % (nome[:9], tr_, es, m) for m in msg]
            k = (nome[6:9], es)
            aperto.setdefault(k, []).append(float(inf / sup))
            acima.setdefault(k, []).append(float(sup / e - 1))
            print("alvo %d eps %s | inferior %s | superior %s | sup/eps-1 %s | inf/sup %s%s"
                  % (tr_, es, mp.nstr(inf, 12), mp.nstr(sup, 12), mp.nstr(sup / e - 1, 3), mp.nstr(inf / sup, 8),
                     "" if not msg else " | FALHA: " + "; ".join(msg)), flush=True)

    print("\n=== RESUMO (50 digitos) ===")
    for k in sorted(aperto):
        a, s = aperto[k], acima[k]
        print("%s eps %s: n = %d | aperto inf/sup mediana %.6f min %.6f max %.6f | superior > eps estrito em %d (max excesso rel %.2e)"
              % (k[0], k[1], len(a), np.median(a), min(a), max(a), sum(x > 0 for x in s), max(s)))
    if falhas:
        print("\nVEREDITO 50 DIGITOS: NAO VERIFICADA (%d falhas)" % len(falhas))
        for m in falhas[:30]:
            print("  ", m)
    else:
        print("\nVEREDITO 50 DIGITOS: VERIFICADA ((a) inferior <= superior e (b) superior <= eps(1+1e-9) em todas)")


if __name__ == "__main__":
    main()
