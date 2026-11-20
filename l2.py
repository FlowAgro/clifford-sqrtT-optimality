"""L2: expoente de denominador de Bloch x contagem de silabas nos circuitos publicados do cyclosynth.

Ver PREREGISTRO_L2.md. Uso: python -X utf8 l2.py
Aritmetica exata em Z[zeta], zeta = e^(2 pi i/16), zeta^8 = -1. Elemento = lista de 8 inteiros.
Matriz U = M / sqrt(2)^h, com M 2x2 sobre Z[zeta] guardada por colunas: M = [[a, c], [b, d]].
"""
import csv
import random
import re
from collections import Counter
from statistics import median

ZERO = [0] * 8
UM = [1] + [0] * 7


def mul(a, b):
    r = [0] * 8
    for i, x in enumerate(a):
        if x:
            for j, y in enumerate(b):
                if y:
                    k = i + j
                    if k >= 8:
                        r[k - 8] -= x * y
                    else:
                        r[k] += x * y
    return r


def add(a, b):
    return [x + y for x, y in zip(a, b)]


def sub(a, b):
    return [x - y for x, y in zip(a, b)]


def neg(a):
    return [-x for x in a]


def zpow(a, k):
    """a * zeta^k."""
    k %= 16
    r = [0] * 8
    for i, x in enumerate(a):
        j = (i + k) % 16
        if j >= 8:
            r[j - 8] -= x
        else:
            r[j] += x
    return r


def conj(a):
    r = [0] * 8
    r[0] = a[0]
    for j in range(1, 8):
        r[8 - j] = -a[j]
    return r


def vP(a):
    """Valoracao no primo (1 - zeta), o unico acima de 2. None para zero."""
    if not any(a):
        return None
    v = 0
    while sum(a) % 2 == 0:
        y0 = (a[0] - sum(a[1:])) // 2
        y = [y0]
        for j in range(1, 8):
            y.append(y[-1] + a[j])
        a = y
        v += 1
    return v


POT = {"Q": 1, "T": 2, "S": 4, "Z": 8, "q": -1, "t": -2, "s": -4}


def matriz(gates):
    """Produto das portas, da esquerda para a direita, como (a, b, c, d, h)."""
    a, b, c, d, h = UM, ZERO, ZERO, UM, 0
    for g in gates:
        if g == "I":  # identidade (src/synthesis/cliffords.rs do cyclosynth); aparece em 3 linhas
            continue
        if g in POT:
            c, d = zpow(c, POT[g]), zpow(d, POT[g])
        elif g == "H":
            a, c = add(a, c), sub(a, c)
            b, d = add(b, d), sub(b, d)
            h += 1
        elif g == "X":
            a, c = c, a
            b, d = d, b
        elif g == "Y":
            a, c = zpow(c, 4), zpow(a, -4)
            b, d = zpow(d, 4), zpow(b, -4)
        else:
            raise ValueError("porta desconhecida: %r" % g)
    return a, b, c, d, h


def mm(A, B):
    (a, b, c, d), (e, f, g, k) = A, B  # colunas: [[a, c], [b, d]]
    return (add(mul(a, e), mul(c, f)), add(mul(b, e), mul(d, f)),
            add(mul(a, g), mul(c, k)), add(mul(b, g), mul(d, k)))


I4 = zpow(UM, 4)
SIG = [(ZERO, UM, UM, ZERO), (ZERO, I4, neg(I4), ZERO), (UM, ZERO, ZERO, neg(UM))]


def expoente_bloch(gates):
    """(N, unitaria): N = maior expoente de denominador em beta = 2cos(pi/8) na matriz de Bloch."""
    a, b, c, d, h = matriz(gates)
    M = (a, b, c, d)
    Md = (conj(a), conj(c), conj(b), conj(d))
    P = mm(M, Md)
    dois_h = [2**h] + [0] * 7
    unit = P[0] == dois_h and P[3] == dois_h and not any(P[1]) and not any(P[2])
    N = 0
    for sj in SIG:
        A11, A21, A12, A22 = mm(mm(M, sj), Md)
        for num in (add(A21, A12), mul(I4, sub(A12, A21)), sub(A11, A22)):
            v = vP(num)
            if v is not None:
                r2 = 8 * (h + 1) - v  # v_P(2) = 8, v_P(beta) = 2
                assert r2 % 2 == 0
                N = max(N, r2 // 2)
    return N, unit


def silabas(gates):
    t = q = 0
    for blk in re.split("[HXY]", gates):
        k = sum(POT.get(ch, 0) for ch in blk) % 4
        if k in (1, 3):
            q += 1
        elif k == 2:
            t += 1
    return t, q


def main():
    # K2
    k2 = {g: expoente_bloch(g)[0] for g in ("T", "Q", "QT", "S", "H")}
    print("K2 portas isoladas:", k2, "->", "OK" if k2 == {"T": 2, "Q": 3, "QT": 3, "S": 0, "H": 0} else "FALHOU")
    linhas = []
    for arq in ("u3", "rz"):
        for r in csv.DictReader(open("dados/%s.csv" % arq, encoding="utf8")):
            linhas.append((arq, r))
    k1 = k4 = 0
    div = Counter()
    exemplos = []
    dados = []
    for arq, r in linhas:
        g = r["gates"]
        t, q = silabas(g)
        if (t, q) != (int(r["t_count"]), int(r["q_count"])) or abs(t + 3 * q - float(r["cost"])) > 1e-9:
            k1 += 1
        N, unit = expoente_bloch(g)
        k4 += not unit
        if N != 2 * t + 3 * q:
            div[(arq, r["synthesizer"], "N<" if N < 2 * t + 3 * q else "N>")] += 1
            if len(exemplos) < 5:
                exemplos.append((arq, r["epsilon"], r["trial"], r["synthesizer"], t, q, N))
        dados.append((arq, r["epsilon"], r["synthesizer"], t, q, N, g))
    print("K1 contador de silabas: divergencias = %d/%d -> %s" % (k1, len(linhas), "OK" if k1 == 0 else "FALHOU"))
    print("K4 unitariedade exata: falhas = %d/%d -> %s" % (k4, len(linhas), "OK" if k4 == 0 else "FALHOU"))
    # K3
    rnd = random.Random(20261114)
    am = rnd.sample([x for x in dados if x[2] == "cyclosynth_sqrt_t"], 200)
    mud = sum(expoente_bloch("HQH" + x[6])[0] != 2 * x[3] + 3 * x[4] for x in am)
    print("K3 mutacao: detectadas = %d/200 -> %s" % (mud, "OK" if mud >= 190 else "FALHOU"))
    # H1
    print("\nH1: linhas com N != 2t + 3q = %d/%d" % (sum(div.values()), len(linhas)))
    for k, v in sorted(div.items()):
        print("   ", k, v)
    for e in exemplos:
        print("    exemplo (arq, eps, alvo, sint, t, q, N):", e)
    # H2 e descritivo
    for arq in ("u3", "rz"):
        print("\n%s.csv, cyclosynth_sqrt_t:" % arq)
        print("  eps    n   mediana q*  max q*  fracao q*<=4  mediana c*/(N*/2)  distribuicao de q*")
        for eps in ("1e-03", "1e-04", "1e-05", "1e-06", "1e-07", "1e-08"):
            xs = [x for x in dados if x[0] == arq and x[1] == eps and x[2] == "cyclosynth_sqrt_t"]
            qs = [x[4] for x in xs]
            raz = [(x[3] + 3 * x[4]) / (x[5] / 2) for x in xs if x[5] > 0]
            print("  %s %d  %5.1f  %5d  %10.3f  %14.3f   %s" % (
                eps, len(xs), median(qs), max(qs), sum(v <= 4 for v in qs) / len(qs), median(raz),
                dict(sorted(Counter(qs).items()))))
        alvo = [x[4] for x in dados if x[0] == arq and x[2] == "cyclosynth_sqrt_t" and x[1] in ("1e-06", "1e-07", "1e-08")]
        print("  eps <= 1e-6: fracao com D <= 8 (q* <= 4) = %d/%d = %.3f" % (
            sum(v <= 4 for v in alvo), len(alvo), sum(v <= 4 for v in alvo) / len(alvo)))


if __name__ == "__main__":
    main()
