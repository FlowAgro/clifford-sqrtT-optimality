"""SONDA A: onde fica o otimo no plano (numero de sqrtT, profundidade)? Um certificado "separado por q" teria chance?

CONTA QUE MOTIVA (Forest et al., Teorema 4.1, conferido na L2): N = 2*n_T + 3*q, logo custo = N/2 + 1,5*q.
O certificado dos autores so fecha quando o horizonte cobre N ~ 2*custo, porque circuitos fundos com POUCAS
sqrtT sao baratos. Se separarmos por q, circuitos com q >= q0 ficam caros mais cedo; os com q < q0 precisam de
outro argumento. O unico solucionador de referencia sem sqrtT e a linha de base Clifford+T (q = 0).

O QUE A SONDA MEDE (busca exaustiva propria da L3, encontro no meio, certificada ate o custo de cobertura):
  C(q) = menor custo de um circuito eps-proximo com EXATAMENTE q silabas sqrtT, para q = 0, 1, 2, ...
  e, para o melhor circuito de cada q, N = 2*n_T + 3*q (Forest) e o lde da matriz (l2.matriz, outra escala).
  Linha de base: custo do cyclosynth com sqrt_t=False (Clifford+T), se a biblioteca estiver disponivel.

REGRAS FIXADAS ANTES DE RODAR (19/11/2026):
  P1 (linha de base serve de certificado para q = 0?): falha se C(0) < custo da linha de base em >= 5% dos alvos.
  P2 (vale separar por q?): seja q_opt o q do otimo. A separacao so ajuda se o otimo tiver q_opt >= 2 em
     >= 50% dos alvos E os circuitos de q pequeno (0, 1) forem mais caros que o otimo por >= 2 unidades
     em >= 50% dos alvos. Se falhar, a ideia A encerra.
  Descritivo: N de Forest e lde da matriz do otimo, por q.

USO (WSL, venv do cyclosynth, dentro de l2_ciclotomico):
  python sonda_A_perfil_q.py A --n 100      # eps = 1e-3, cobertura custo <= 21 (~2 min de tabela + ~1 s/alvo)
  python sonda_A_perfil_q.py B --n 40       # eps = 1e-4, cobertura custo <= 27 (tabela maior; ~10-20 min)
Log automatico em logs/sonda_A_<modo>_<data>.log e tabela em logs/sonda_A_<modo>_<data>.csv.
"""
import argparse
import cmath
import csv
import os
import time
from collections import Counter

import numpy as np

import l2
import l3
import registro

CONF = {"A": dict(eps="1e-03", a=12, b=11), "B": dict(eps="1e-04", a=16, b=13)}


def preparar(a, b, eps):
    al, be, par, sil, custo = l3.gerar(a)
    nR = int(np.searchsorted(custo, b, side="right"))
    cl = l3.cliffords()
    assert len(cl) == 24, "esperava 24 Cliffords"
    qR = np.empty((nR * 24, 4))
    for i, s in enumerate(cl):
        ca, cb = l3.su2(s)
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
    return dict(al=al, be=be, par=par, sil=sil, custo=custo, nR=nR, cl=cl, qR=qR, cR=cR, r=r, hc=hc,
                pesos=pesos, KR=chave[ordem], ordem=ordem, eps=eps)


def candidatos(T, V, teto):
    """Todos os pares (custo_total, i_esq, i_dir) com custo_total <= teto e quase-proximidade (mesmo laco da L3)."""
    al, be, KR, ordem, qR, cR, custo, pesos = T["al"], T["be"], T["KR"], T["ordem"], T["qR"], T["cR"], T["custo"], T["pesos"]
    r, hc, eps = T["r"], T["hc"], T["eps"]
    Vs = V / cmath.sqrt(np.linalg.det(V))
    va, vb = Vs[0, 0], Vs[1, 0]
    out = []
    BL = 1_000_000
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
                    rows, cel = np.arange(len(pp)), lo
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
                ok = L < len(KR)
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
                    if ctot <= teto:
                        out.append((ctot, s0 + int(rr[j]), int(ri[j])))
    return sorted(set(out))


def perfil(T, V, cob):
    """{q: (custo, N_Forest, lde_matriz, gates)} com o menor custo reconferido em aritmetica exata para cada q."""
    best = {}
    for ctot, li, ri in candidatos(T, V, cob):
        g = l3.palavra(li, T["par"], T["sil"]) + l3.palavra(ri % T["nR"], T["par"], T["sil"]) + T["cl"][ri // T["nR"]]
        t, q = l2.silabas(g)
        c = t + 3 * q
        if q in best and best[q][0] <= c:
            continue
        if l3.distancia_exata(g, V) < T["eps"]:
            best[q] = (c, 2 * t + 3 * q, l2.matriz(g)[4], g)
    return best


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("modo", choices=["A", "B"])
    p.add_argument("--n", type=int, default=100)
    p.add_argument("--primeiro", type=int, default=300, help="primeiro trial (padrao 300: fora dos usados na L3 em B)")
    args = p.parse_args()
    log = registro.iniciar_log("sonda_A_" + args.modo)
    cf = CONF[args.modo]
    eps = float(cf["eps"])
    cob = cf["a"] + cf["b"] - 2
    try:
        import cyclosynth
        base = cyclosynth.Synthesizer(epsilon=eps, sqrt_t=False)
        print("linha de base Clifford+T: cyclosynth sqrt_t=False")
    except Exception as e:  # sem a biblioteca, a P1 fica sem teste
        base = None
        print("AVISO: cyclosynth indisponivel (%s); P1 nao sera testada." % e)

    t0 = time.time()
    T = preparar(cf["a"], cf["b"], eps)
    print("tabela pronta em %.0f s: E = %d palavras, D x Clifford = %d pontos; cobertura certificada: custo <= %d"
          % (time.time() - t0, len(T["al"]), len(T["KR"]), cob), flush=True)

    alvos = {}
    for row in csv.DictReader(open(os.path.join(l3.__file__.rsplit(os.sep, 1)[0], "dados", "u3.csv"), encoding="utf8")):
        tr = int(row["trial"])
        if row["epsilon"] == cf["eps"] and row["synthesizer"] == "cyclosynth_sqrt_t" and args.primeiro <= tr < args.primeiro + args.n:
            alvos[tr] = row
    csv_out = log[:-4] + ".csv"
    w = csv.writer(open(csv_out, "w", encoding="utf8", newline=""))
    w.writerow(["trial", "custo_cyclosynth", "custo_base_CT", "opt", "q_opt", "N_Forest_opt", "lde_opt", "perfil_q:custo:N:lde"])

    res = []
    for tr in sorted(alvos):
        row = alvos[tr]
        V = l3.alvo(float(row["alpha"]), float(row["beta"]), float(row["gamma"]))
        pf = perfil(T, V, cob)
        cb = None
        if base is not None:
            rb = base.synthesize_u3(float(row["beta"]), float(row["alpha"]), float(row["gamma"]))
            cb = int(rb.cost) if rb else None
        if not pf:
            print("alvo %d: nada com custo <= %d (fora da cobertura)" % (tr, cob))
            continue
        opt = min(v[0] for v in pf.values())
        q_opt = min(q for q, v in pf.items() if v[0] == opt)
        N_opt, lde_opt = pf[q_opt][1], pf[q_opt][2]
        txt = " ".join("q%d:%d:N%d:lde%d" % (q, v[0], v[1], v[2]) for q, v in sorted(pf.items()))
        print("alvo %d: cyclosynth %s | base C+T %s | OPT %d (q=%d, N=%d, lde=%d) | perfil %s"
              % (tr, row["cost"], cb, opt, q_opt, N_opt, lde_opt, txt), flush=True)
        w.writerow([tr, row["cost"], cb, opt, q_opt, N_opt, lde_opt, txt])
        res.append((tr, int(float(row["cost"])), cb, opt, q_opt, N_opt, pf, lde_opt))

    n = len(res)
    print("\n=== RESUMO (eps = %s, n = %d alvos, %.0f s) ===" % (cf["eps"], n, time.time() - t0))
    if not n:
        return
    print("Distribuicao de q no otimo:", dict(sorted(Counter(r[4] for r in res).items())))
    print("Otimo: N de Forest medio %.1f | lde da matriz medio %.1f | OPT medio %.2f | cyclosynth medio %.2f"
          % (np.mean([r[5] for r in res]), np.mean([r[7] for r in res]), np.mean([r[3] for r in res]),
             np.mean([r[1] for r in res])))
    # P1
    if base is not None:
        com0 = [r for r in res if 0 in r[6] and r[2] is not None]
        piores = sum(r[6][0][0] < r[2] for r in com0)
        print("P1: alvos com circuito q=0 na cobertura: %d | base C+T acima do melhor q=0: %d (%.1f%%) -> %s"
              % (len(com0), piores, 100 * piores / max(1, len(com0)),
                 "NAO TESTAVEL (nenhum q=0 na cobertura)" if not com0 else
                 ("FALHOU" if piores >= 0.05 * len(com0) else "passou")))
        if com0:
            print("    diferenca media base - C(0): %.2f" % np.mean([r[2] - r[6][0][0] for r in com0]))
    # P2
    q2 = sum(r[4] >= 2 for r in res)
    gap = []
    for r in res:
        baixos = [r[6][q][0] for q in (0, 1) if q in r[6]]
        gap.append((min(baixos) - r[3]) if baixos else None)  # None = q<=1 fora da cobertura (mais caro que cob)
    sep = sum(1 for g in gap if g is None or g >= 2)
    ok2 = q2 >= 0.5 * n and sep >= 0.5 * n
    print("P2: otimo com q >= 2 em %d/%d (%.0f%%); q in {0,1} mais caro que o otimo por >= 2 em %d/%d (%.0f%%) -> %s"
          % (q2, n, 100 * q2 / n, sep, n, 100 * sep / n, "passou" if ok2 else "FALHOU"))
    print("    (None conta como 'mais caro': nenhum circuito com q <= 1 dentro da cobertura)")
    print("VEREDITO A:", "seguir para o desenho do certificado" if ok2 else "a separacao por q nao ajuda: encerrar A")
    print("tabela por alvo em:", csv_out)


if __name__ == "__main__":
    main()
