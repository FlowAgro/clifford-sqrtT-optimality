"""L3: busca exaustiva por encontro no meio de circuitos Clifford+sqrt(T) mais baratos que os do cyclosynth.

Ver PREREGISTRO_L3.md. Uso: python -X utf8 l3.py A | B  [primeiro ultimo]
Unitarias em SU(2) como pares (alfa, beta): [[alfa, -conj(beta)], [beta, conj(alfa)]].
Todo circuito achado e reconferido em aritmetica exata com l2.matriz e custeado com l2.silabas.
"""
import cmath
import csv
import math
import sys
import time

import numpy as np

import l2

CONF = {"A": dict(eps="1e-03", a=12, b=11, alvos=(0, 499)),
        "B": dict(eps="1e-04", a=16, b=13, alvos=(0, 149))}

R2 = math.sqrt(0.5)
PORTA = {
    "I": np.eye(2, dtype=complex),
    "H": np.array([[R2, R2], [R2, -R2]], dtype=complex),
    "X": np.array([[0, 1], [1, 0]], dtype=complex),
    "Y": np.array([[0, -1j], [1j, 0]], dtype=complex),
}
for nome, k in l2.POT.items():
    PORTA[nome] = np.array([[1, 0], [0, cmath.exp(1j * math.pi * k / 8)]], dtype=complex)


def su2(gates):
    M = np.eye(2, dtype=complex)
    for g in gates:
        M = M @ PORTA[g]
    M = M / cmath.sqrt(np.linalg.det(M))
    return M[0, 0], M[1, 0]


ZS = {1: "Q", 2: "T", 3: "QT"}
SIL_TXT = {}
SIL = {}
for p in range(3):
    for a in (1, 2, 3):
        txt = ("H%sH", "SH%sHs", "%s")[p] % ZS[a]
        SIL_TXT[p * 3 + a - 1] = txt
        SIL[(p, a)] = su2(txt)


def cliffords():
    vistos, fila = {}, [""]
    while fila:
        s = fila.pop(0)
        al, be = su2(s)
        v = np.array([al.real, al.imag, be.real, be.imag])
        i = int(np.argmax(np.abs(v) > 1e-9))
        v = v * (1 if v[i] > 0 else -1)
        ch = tuple(np.round(v, 6) + 0.0)
        if ch in vistos:
            continue
        vistos[ch] = s
        fila += [s + "H", s + "S"]
    return sorted(vistos.values(), key=lambda s: (len(s), s))


def gerar(amax):
    """Todas as palavras canonicas de custo <= amax, em ordem de custo."""
    niv = {0: (np.array([1 + 0j]), np.array([0j]), np.array([-1], dtype=np.int8),
               np.array([-1], dtype=np.int64), np.array([-1], dtype=np.int8))}
    off = {0: 0}
    total = 1
    for j in range(1, amax + 1):
        partes = []
        for src, aops in ((j - 1, (2,)), (j - 3, (1, 3))):
            if src < 0:
                continue
            al, be, last, _, _ = niv[src]
            for p in range(3):
                m = np.nonzero(last != p)[0]
                for a in aops:
                    sa, sb = SIL[(p, a)]
                    partes.append((al[m] * sa - np.conj(be[m]) * sb, be[m] * sa + np.conj(al[m]) * sb,
                                   np.full(len(m), p, dtype=np.int8), off[src] + m,
                                   np.full(len(m), p * 3 + a - 1, dtype=np.int8)))
        niv[j] = tuple(np.concatenate([x[i] for x in partes]) for i in range(5))
        off[j] = total
        total += len(niv[j][0])
    al, be, _, par, sil = (np.concatenate([niv[j][i] for j in range(amax + 1)]) for i in range(5))
    custo = np.concatenate([np.full(len(niv[j][0]), j, dtype=np.int8) for j in range(amax + 1)])
    return al, be, par, sil, custo


def palavra(i, par, sil):
    out = []
    while par[i] >= 0:
        out.append(SIL_TXT[int(sil[i])])
        i = par[i]
    return "".join(reversed(out))


def alvo(alpha, beta, gamma):
    def rz(t):
        return np.array([[cmath.exp(-0.5j * t), 0], [0, cmath.exp(0.5j * t)]])
    ry = np.array([[math.cos(beta / 2), -math.sin(beta / 2)], [math.sin(beta / 2), math.cos(beta / 2)]], dtype=complex)
    return rz(alpha) @ ry @ rz(gamma)


def distancia_exata(gates, V):
    a, b, c, d, h = l2.matriz(gates)
    z = [cmath.exp(1j * math.pi * k / 8) for k in range(8)]
    ev = lambda e: sum(x * w for x, w in zip(e, z))
    tr = (ev(a).conjugate() * V[0, 0] + ev(b).conjugate() * V[1, 0]
          + ev(c).conjugate() * V[0, 1] + ev(d).conjugate() * V[1, 1]) / 2 ** (h / 2)
    return math.sqrt(max(0.0, 1 - abs(tr) ** 2 / 4))


def main():
    modo = sys.argv[1]
    cf = CONF[modo]
    eps = float(cf["eps"])
    a, b = cf["a"], cf["b"]
    ini, fim = (int(sys.argv[2]), int(sys.argv[3])) if len(sys.argv) > 3 else cf["alvos"]
    t0 = time.time()
    al, be, par, sil, custo = gerar(a)
    nR = int(np.searchsorted(custo, b, side="right"))
    cum13 = int(np.searchsorted(custo, min(13, a), side="right"))
    print("T3 contagem: palavras com custo <= %d: %d (esperado 403660 para 13; 155596 para 12)" % (min(13, a), cum13))
    print("palavras E (custo <= %d): %d ; palavras D (custo <= %d): %d" % (a, len(al), b, nR), flush=True)
    cl = cliffords()
    print("T4 Cliffords: %d distintos; custo maximo pela regra dos autores: %d -> %s"
          % (len(cl), max(sum(l2.silabas(s)) for s in cl), "OK" if len(cl) == 24 and all(sum(l2.silabas(s)) == 0 for s in cl) else "FALHOU"))
    # tabela da direita: D . Clifford
    qR = np.empty((nR * 24, 4))
    for i, s in enumerate(cl):
        ca, cb = su2(s)
        ra = al[:nR] * ca - np.conj(be[:nR]) * cb
        rb = be[:nR] * ca + np.conj(al[:nR]) * cb
        blk = qR[i * nR:(i + 1) * nR]
        blk[:, 0], blk[:, 1], blk[:, 2], blk[:, 3] = ra.real, ra.imag, rb.real, rb.imag
    cR = np.tile(custo[:nR], 24)
    r = 1.001 * eps
    hc = 40 * r
    Mb = int(2 / hc) + 3
    pesos = np.array([Mb**3, Mb**2, Mb, 1], dtype=np.int64)
    chave = np.floor((qR + 1) / hc).astype(np.int64) @ pesos
    ordem = np.argsort(chave, kind="stable")
    KR = chave[ordem]
    del chave
    print("tabela D pronta: %d pontos, %.0f s" % (len(KR), time.time() - t0), flush=True)

    linhas = {}
    for row in csv.DictReader(open("dados/u3.csv", encoding="utf8")):
        if row["epsilon"] == cf["eps"] and row["synthesizer"] == "cyclosynth_sqrt_t" and ini <= int(row["trial"]) <= fim:
            linhas[int(row["trial"])] = row
    cob = a + b - 2
    cont = {"SUBOTIMO": 0, "CERTIFICADO": 0, "INDETERMINADO": 0}
    t1_falhas = t2_falhas = 0
    economia = []
    saida = open("l3_%s_testemunhas.csv" % modo, "a", encoding="utf8", newline="")
    wr = csv.writer(saida)
    BL = 1_000_000
    for trial in sorted(linhas):
        row = linhas[trial]
        V = alvo(float(row["alpha"]), float(row["beta"]), float(row["gamma"]))
        cstar = int(float(row["cost"]))
        d_aut = distancia_exata(row["gates"], V)
        if abs(d_aut - float(row["distance"])) > 1e-4 * float(row["distance"]):
            t1_falhas += 1
        Vs = V / cmath.sqrt(np.linalg.det(V))
        va, vb = Vs[0, 0], Vs[1, 0]
        achados = []
        for s0 in range(0, len(al), BL):
            la, lb = al[s0:s0 + BL], be[s0:s0 + BL]
            pa = np.conj(la) * va + np.conj(lb) * vb
            pb = -lb * va + la * vb
            p = np.stack([pa.real, pa.imag, pb.real, pb.imag], axis=1)
            for sg in (1.0, -1.0):
                pp = p * sg
                lo = np.floor((pp - r + 1) / hc).astype(np.int64)
                hi = np.floor((pp + r + 1) / hc).astype(np.int64)
                dif = hi != lo
                for combo in range(16):
                    bits = [(combo >> i) & 1 for i in range(4)]
                    if combo == 0:
                        rows = np.arange(len(pp))
                        cel = lo
                    else:
                        m = np.ones(len(pp), dtype=bool)
                        for i in range(4):
                            if bits[i]:
                                m &= dif[:, i]
                        rows = np.nonzero(m)[0]
                        if len(rows) == 0:
                            continue
                        cel = np.where(np.array(bits, dtype=bool)[None, :], hi[rows], lo[rows])
                    key = cel @ pesos
                    L = np.searchsorted(KR, key)
                    ok = (L < len(KR))
                    ok[ok] = KR[L[ok]] == key[ok]
                    if not ok.any():
                        continue
                    rows, key, L = rows[ok], key[ok], L[ok]
                    Rr = np.searchsorted(KR, key, side="right")
                    n = Rr - L
                    rr = np.repeat(rows, n)
                    pos = np.repeat(L, n) + (np.arange(n.sum()) - np.repeat(np.cumsum(n) - n, n))
                    ri = ordem[pos]
                    dot = np.einsum("ij,ij->i", pp[rr], qR[ri])
                    h = np.nonzero(1 - dot * dot < eps * eps * (1 + 1e-6))[0]
                    for j in h:
                        ctot = int(custo[s0 + rr[j]]) + int(cR[ri[j]])
                        if ctot <= cstar:
                            achados.append((ctot, s0 + int(rr[j]), int(ri[j])))
        melhor = None
        for ctot, li, ri in sorted(set(achados)):
            g = palavra(li, par, sil) + palavra(ri % nR, par, sil) + cl[ri // nR]
            d = distancia_exata(g, V)
            if d < eps:
                cv = sum(x * w for x, w in zip(l2.silabas(g), (1, 3)))
                if melhor is None or (cv, d) < melhor[:2]:
                    melhor = (cv, d, g)
        m = melhor[0] if melhor else None
        if cstar <= cob and (m is None or m > cstar):
            t2_falhas += 1
        if m is not None and m < cstar:
            ver = "SUBOTIMO"
            economia.append(cstar - m)
            wr.writerow([cf["eps"], trial, cstar, m, "%.6e" % melhor[1], melhor[2], row["gates"]])
            saida.flush()
        elif cstar - 1 <= cob:
            ver = "CERTIFICADO"
        else:
            ver = "INDETERMINADO"
        cont[ver] += 1
        print("alvo %d: c* = %d, melhor achado = %s, candidatos = %d -> %s%s" % (
            trial, cstar, m, len(set(achados)), ver,
            (" | d = %.4e | %s" % (melhor[1], melhor[2])) if ver == "SUBOTIMO" else ""), flush=True)
    n = len(linhas)
    print("\nconfiguracao %s: eps = %s, a = %d, b = %d, cobertura = custo <= %d, alvos %d a %d (n = %d), %.0f s"
          % (modo, cf["eps"], a, b, cob, ini, fim, n, time.time() - t0))
    print("T1 convencao (distancia dos autores reproduzida): falhas = %d/%d -> %s" % (t1_falhas, n, "OK" if t1_falhas == 0 else "FALHOU"))
    print("T2 cobertura (acha custo <= c* quando c* <= %d): falhas = %d -> %s" % (cob, t2_falhas, "OK" if t2_falhas == 0 else "FALHOU"))
    print("veredito:", cont)
    if economia:
        print("economia nos subotimos: media = %.2f, max = %d (em estados T)" % (sum(economia) / len(economia), max(economia)))


if __name__ == "__main__":
    main()
