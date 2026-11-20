"""Conferencia independente das testemunhas da L3.

Custo: funcao dos PROPRIOS autores (dados/_cost_autores.py = scripts/_cost.py do cyclosynth, sem alteracao).
Distancia: mpmath com 60 digitos, matrizes de porta em ponto flutuante de alta precisao (outro caminho que
nao passa por l2.matriz). Uso: python -X utf8 l3_conferir.py A B
"""
import csv
import sys

import mpmath as mp

sys.path.insert(0, "dados")
import _cost_autores as aut

mp.mp.dps = 60
r2 = mp.sqrt(mp.mpf(1) / 2)
G = {"I": mp.eye(2), "H": mp.matrix([[r2, r2], [r2, -r2]]), "X": mp.matrix([[0, 1], [1, 0]]),
     "Y": mp.matrix([[0, -1j], [1j, 0]])}
for nome, k in {"Q": 1, "T": 2, "S": 4, "Z": 8, "q": -1, "t": -2, "s": -4}.items():
    G[nome] = mp.matrix([[1, 0], [0, mp.expjpi(mp.mpf(k) / 8)]])


def unit(gates):
    M = mp.eye(2)
    for g in gates:
        M = M * G[g]
    return M


def alvo(a, b, g):
    rz = lambda t: mp.matrix([[mp.expj(-t / 2), 0], [0, mp.expj(t / 2)]])
    ry = mp.matrix([[mp.cos(b / 2), -mp.sin(b / 2)], [mp.sin(b / 2), mp.cos(b / 2)]])
    return rz(a) * ry * rz(g)


def dist(gates, V):
    U = unit(gates)
    tr = sum(mp.conj(U[i, j]) * V[i, j] for i in range(2) for j in range(2))
    return mp.sqrt(1 - abs(tr) ** 2 / 4)


ang = {}
for r in csv.DictReader(open("dados/u3.csv", encoding="utf8")):
    ang[(r["epsilon"], r["trial"])] = tuple(mp.mpf(r[k]) for k in ("alpha", "beta", "gamma"))
for modo in sys.argv[1:]:
    n = ok = 0
    pior = 0
    for eps, trial, cstar, m, d, g, g_aut in csv.reader(open("l3_%s_testemunhas.csv" % modo, encoding="utf8")):
        V = alvo(*ang[(eps, trial)])
        c_nosso, c_deles = aut.syllable_cost(g), aut.syllable_cost(g_aut)
        d_nosso, d_deles = dist(g, V), dist(g_aut, V)
        e = mp.mpf(eps)
        bom = c_nosso == int(m) and c_deles == int(cstar) and c_nosso < c_deles and d_nosso < e and d_deles < e
        n += 1
        ok += bom
        pior = max(pior, d_nosso / e)
        if not bom:
            print("FALHA", modo, trial, c_nosso, c_deles, mp.nstr(d_nosso, 8), mp.nstr(d_deles, 8))
    print("%s: testemunhas conferidas = %d/%d (custo pela funcao dos autores, distancia com 60 digitos); maior d/eps = %s"
          % (modo, ok, n, mp.nstr(pior, 6)))
