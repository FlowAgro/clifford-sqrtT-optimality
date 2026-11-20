"""VERIFICACAO INDEPENDENTE da sonda M (mistura Clifford+sqrtT), passo 2 de RETOMAR_MISTURA.md.

Nao importa l2, l3 nem sonda_M_mistura: tudo e refeito do zero a partir das cadeias de portas e das
probabilidades gravadas nos CSVs (as p estao arredondadas a 1e-6; verificamos a mistura ARREDONDADA, que e a
que alguem usaria, renormalizada para somar 1).

O que e conferido, por mistura (alvo x eps efetivo):
  (i)   ARITMETICA EXATA: cada circuito e multiplicado em Z[zeta16] (inteiros, zeta^8 = -1), porta a porta, com
        matrizes 2x2 genericas (nao o atalho por colunas do l2). Confere U U^dag = 2^h I EXATAMENTE no anel.
  (ii)  DISTANCIA: d = sqrt(1 - |tr(U^dag V)|^2 / 4) contra o alvo V = Rz(alpha) Ry(beta) Rz(gamma).
        Controle de convencao (alvo e ordem das portas): o mesmo calculo no circuito PUBLICADO do
        cyclosynth_sqrt_t tem de reproduzir a coluna 'distance' de dados/u3.csv (tolerancia relativa 1e-3,
        a coluna tem 7 algarismos).
  (iii) CUSTO: recalculado por dados/_cost_autores.py (espelho do gates_cost em Rust dos autores); o custo
        esperado sum p_i custo_i tem de bater com o CSV (tolerancia 1e-3).
  (iv)  COTA SUPERIOR recalculada: |sum p_i c_i a_i| + sum p_i |a_i|^2 <= eps_ef (1 + 1e-6).
  (v)   COTA INFERIOR da distancia diamante (1/2)||E - Id||_diamante, por outro caminho: maximo, sobre estados
        puros psi de qubit + ancila, de (1/2)||(E - Id) x I (psi psi^dag)||_1 (cada psi da uma cota inferior
        valida; busca com 24 partidas aleatorias + o estado maximamente emaranhado, Nelder-Mead).
        Tem de valer: cota inferior <= eps_ef (senao a cota superior esta ERRADA).
  CONTROLE do otimizador: para um circuito isolado (canal unitario W = c I - i a.sigma) a distancia diamante
        e exatamente |a|. A busca, aplicada ao circuito de maior |a| de cada suporte, tem de achar >= 0,99 |a|.

REGRA FIXADA ANTES DE RODAR (19/11/2026): a sonda M fica VERIFICADA se (i)-(v) e o controle valem em TODAS
as misturas viaveis dos dois CSVs oficiais (sonda_M_A_20261119_174254, sonda_M_B_20261119_174751) e o
controle de convencao vale em todos os alvos. Qualquer falha = nao verificada, e o motivo e impresso.
Descritivo: aperto = cota inferior / cota superior (perto de 1: a cota e justa; perto de 0: folga escondida).

AJUSTE (19/11, 19h20, ANTES da rodada valida): a 1a tentativa (log verifica_M_20261119_191234, parada em 12
misturas) falhou em (iv)/(v) por 0,1-0,2%: os CSVs oficiais gravaram p com 6 casas e a mistura otima fica
EXATAMENTE em cima de eps, entao o arredondamento a empurra para fora. Correcao so no DADO: a sonda M passou a
gravar p com precisao total (%r) e foi rodada de novo nos mesmos alvos (310+, mesmos modos). Regra inalterada.
A propria 1a tentativa ja mostrou aperto 0,77-1,0: a cota superior e quase justa.

USO (WSL, venv do cyclosynth; precisa so de numpy e scipy):
  python verifica_M.py [csv_A csv_B]   # nomes em logs/; sem argumento, os dois CSVs oficiais
"""
import csv
import math
import os
import sys

import numpy as np
from scipy.optimize import minimize

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(AQUI, "dados"))
from _cost_autores import syllable_cost  # noqa: E402

import registro  # noqa: E402  (so o espelho da tela para o log)

CSVS = ["sonda_M_A_20261119_174254.csv", "sonda_M_B_20261119_174751.csv"]

# ---------- Z[zeta16]: listas de 8 inteiros, zeta^8 = -1 ----------
def r_mul(x, y):
    z = [0] * 8
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                if b:
                    k = i + j
                    if k < 8:
                        z[k] += a * b
                    else:
                        z[k - 8] -= a * b
    return z


def r_add(x, y):
    return [a + b for a, b in zip(x, y)]


def r_conj(x):  # zeta^k -> zeta^-k = -zeta^(8-k)
    z = [x[0]] + [0] * 7
    for k in range(1, 8):
        z[8 - k] -= x[k]
    return z


def r_zeta(k):
    k %= 16
    z = [0] * 8
    if k < 8:
        z[k] = 1
    else:
        z[k - 8] = -1
    return z


ZERO, ONE = [0] * 8, r_zeta(0)
NEG = r_zeta(8)
POT = {"Q": 1, "T": 2, "S": 4, "Z": 8, "q": -1, "t": -2, "s": -4}
GATE = {g: ([[ONE, ZERO], [ZERO, r_zeta(k)]], 0) for g, k in POT.items()}
GATE["H"] = ([[ONE, ONE], [ONE, NEG]], 1)          # (1/sqrt2) [[1,1],[1,-1]]
GATE["X"] = ([[ZERO, ONE], [ONE, ZERO]], 0)
GATE["Y"] = ([[ZERO, r_zeta(12)], [r_zeta(4), ZERO]], 0)  # [[0,-i],[i,0]]
GATE["I"] = ([[ONE, ZERO], [ZERO, ONE]], 0)


def m_mul(A, B):
    return [[r_add(r_mul(A[i][0], B[0][j]), r_mul(A[i][1], B[1][j])) for j in range(2)] for i in range(2)]


def exata(gates):
    """Produto G1 G2 ... Gn (esquerda para a direita) e o numero h de H (fator 2^(-h/2))."""
    M, h = [[ONE, ZERO], [ZERO, ONE]], 0
    for g in gates:
        G, dh = GATE[g]
        M, h = m_mul(M, G), h + dh
    return M, h


def unitaria_exata(M, h):
    Md = [[r_conj(M[j][i]) for j in range(2)] for i in range(2)]
    P = m_mul(M, Md)
    dois_h = [2 ** h] + [0] * 7
    return P[0][0] == dois_h and P[1][1] == dois_h and P[0][1] == ZERO and P[1][0] == ZERO


Z16 = np.exp(1j * np.pi * np.arange(8) / 8)


def numerica(M, h):
    return np.array([[np.dot(M[i][j], Z16) for j in range(2)] for i in range(2)]) / 2 ** (h / 2)


# ---------- alvo, distancia, canal ----------
def alvo(al, be, ga):
    rz = lambda t: np.diag([np.exp(-0.5j * t), np.exp(0.5j * t)])
    ry = np.array([[math.cos(be / 2), -math.sin(be / 2)], [math.sin(be / 2), math.cos(be / 2)]], dtype=complex)
    return rz(al) @ ry @ rz(ga)


def dist(U, V):
    return math.sqrt(max(0.0, 1 - abs(np.trace(U.conj().T @ V)) ** 2 / 4))


SIG = [np.array([[0, 1], [1, 0]], complex), np.array([[0, -1j], [1j, 0]], complex), np.diag([1, -1]).astype(complex)]


def ca(U, V):
    W = V.conj().T @ U
    W = W / np.sqrt(np.linalg.det(W))
    c = (np.trace(W) / 2).real
    a = np.array([(0.5j * np.trace(W @ P)).real for P in SIG])
    return (c, a) if c >= 0 else (-c, -a)


def kraus(c, a):
    return c * np.eye(2) - 1j * sum(x * P for x, P in zip(a, SIG))


def f_diamante(Ws, p, x):
    """(1/2)||(E - Id) x I (psi psi^dag)||_1 com psi dado por 8 reais."""
    psi = (x[:4] + 1j * x[4:]).reshape(2, 2)  # psi[sistema, ancila]
    n = np.linalg.norm(psi)
    if n < 1e-12:
        return 0.0
    psi = psi / n
    rho = psi.reshape(4, 1) @ psi.reshape(1, 4).conj()
    D = -rho
    for W, pi in zip(Ws, p):
        K = np.kron(W, np.eye(2))
        D = D + pi * (K @ rho @ K.conj().T)
    return 0.5 * np.abs(np.linalg.eigvalsh((D + D.conj().T) / 2)).sum()


def inferior(Ws, p, rng, partidas=24):
    bell = np.array([1, 0, 0, 1, 0, 0, 0, 0], float)
    xs = [bell] + [rng.normal(size=8) for _ in range(partidas)]
    melhor = 0.0
    for x0 in xs:
        r = minimize(lambda x: -f_diamante(Ws, p, x), x0, method="Nelder-Mead",
                     options=dict(maxiter=4000, xatol=1e-10, fatol=1e-18))
        melhor = max(melhor, -r.fun, f_diamante(Ws, p, x0))
    return melhor


def main():
    registro.iniciar_log("verifica_M")
    rng = np.random.default_rng(20261119)
    u3 = list(csv.DictReader(open(os.path.join(AQUI, "dados", "u3.csv"), encoding="utf8")))
    pub = {(r["epsilon"], int(r["trial"])): r for r in u3 if r["synthesizer"] == "cyclosynth_sqrt_t"}
    falhas, apertos, n_mix, conv_ok = [], [], 0, 0
    eps_pool = {"sonda_M_A": "1e-03", "sonda_M_B": "1e-04"}
    for nome in (sys.argv[1:] or CSVS):
        pool_eps = eps_pool[nome[:9]]
        linhas = list(csv.DictReader(open(os.path.join(AQUI, "logs", nome), encoding="utf8")))
        print("\n== %s: %d linhas ==" % (nome, len(linhas)))
        vistos = set()
        for r in linhas:
            tr, e = int(r["trial"]), float(r["eps_ef"])
            ref = pub[(pool_eps, tr)]
            V = alvo(float(ref["alpha"]), float(ref["beta"]), float(ref["gamma"]))
            if tr not in vistos:  # controle de convencao no circuito publicado
                vistos.add(tr)
                M, h = exata(ref["gates"])
                dp, dref = dist(numerica(M, h), V), float(ref["distance"])
                if abs(dp - dref) <= 1e-3 * dref and syllable_cost(ref["gates"]) == float(ref["cost"]):
                    conv_ok += 1
                else:
                    falhas.append("convencao alvo %d: d %.6e vs %.6e" % (tr, dp, dref))
            if not r["custo_esperado"]:
                continue
            n_mix += 1
            itens = [s.rsplit(":", 1) for s in r["circuitos_e_probabilidades"].split()]
            gs, p = [g for g, _ in itens], np.array([float(x) for _, x in itens])
            p = p / p.sum()
            Ws, cs, a2, custos, amax = [], np.zeros(3), 0.0, [], 0.0
            for g, pi in zip(gs, p):
                M, h = exata(g)
                if not unitaria_exata(M, h):
                    falhas.append("alvo %d eps %g: U U^dag != 2^h I (%s)" % (tr, e, g))
                U = numerica(M, h)
                c, a = ca(U, V)
                if abs(np.linalg.norm(a) - dist(U, V)) > 1e-12:
                    falhas.append("alvo %d: |a| != d" % tr)
                Ws.append(kraus(c, a))
                cs += pi * c * a
                a2 += pi * (a @ a)
                custos.append(syllable_cost(g))
                if np.linalg.norm(a) > amax:
                    amax, Wmax = float(np.linalg.norm(a)), Ws[-1]
            sup = float(np.linalg.norm(cs) + a2)
            cesp = float(p @ np.array(custos))
            inf = inferior(Ws, p, rng)
            ctrl = inferior([Wmax], [1.0], rng, partidas=4)
            ok = []
            if abs(cesp - float(r["custo_esperado"])) > 1e-3:
                ok.append("custo %.4f vs %s" % (cesp, r["custo_esperado"]))
            if sup > e * (1 + 1e-6):
                ok.append("cota superior %.3e > eps" % sup)
            if inf > e:
                ok.append("COTA INFERIOR %.3e > eps" % inf)
            if ctrl < 0.99 * amax:
                ok.append("controle: achou %.3e de |a| = %.3e" % (ctrl, amax))
            falhas += ["alvo %d eps %g: %s" % (tr, e, s) for s in ok]
            apertos.append((nome[7:9], e, inf / sup))
            print("alvo %d eps %g: custo %.3f | diamante em [%.3e, %.3e] | aperto %.3f | controle %.4f%s"
                  % (tr, e, cesp, inf, sup, inf / sup, ctrl / amax, "" if not ok else " | FALHA: " + "; ".join(ok)),
                  flush=True)

    print("\n=== RESUMO ===")
    print("controle de convencao (circuito publicado reproduz 'distance' e 'cost'): %d alvos ok" % conv_ok)
    print("misturas verificadas: %d" % n_mix)
    for m in ("_A", "_B"):
        for e in sorted({x[1] for x in apertos if x[0] == m}):
            v = [x[2] for x in apertos if x[0] == m and x[1] == e]
            print("modo %s eps %g: aperto (inferior/superior) mediana %.3f, min %.3f, max %.3f (n = %d)"
                  % (m[1], e, np.median(v), min(v), max(v), len(v)))
    if falhas:
        print("\nVEREDITO: NAO VERIFICADA (%d falhas)" % len(falhas))
        for s in falhas[:30]:
            print("  ", s)
    else:
        print("\nVEREDITO: VERIFICADA (todas as condicoes (i)-(v), controle do otimizador e de convencao)")


if __name__ == "__main__":
    main()
