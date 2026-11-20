"""L5 (EXPLORATORIO, sem pre-registro): lei de contagem para o custo otimo de sintese de 1 qubit.

Para Clifford + rotacoes z do nivel l da hierarquia de Clifford (l = 3 e a porta T), com preco c_l por
silaba e 2^(l-3) classes de silaba por nivel (forma canonica unica: Forest et al., arXiv:1501.04944),
o numero de unitarias com custo <= C cresce como x^C, x raiz de  1 = 2 * soma_l 2^(l-3) * x^(-c_l).
Limite inferior de custo para alvo tipico: (3 / log2 x) * log2(1/eps). Sem parametro livre.

Comparacao com os ajustes MEDIDOS por busca exaustiva em Mooney, Hill, Hollenberg, Quantum 5, 396 (2021),
arXiv:2005.05581, Tabelas 1 e 3 (inclinacao por log10(1/eps)). Antes de calcular eu ja tinha visto as
reducoes percentuais deles (34/43/49/54 e 29/31/31/31), nao os precos nem as inclinacoes.
"""
import math


def raiz(precos):
    f = lambda x: 2 * sum(m * x ** (-c) for m, c in precos) - 1
    lo, hi = 1.0000001, 1e6
    for _ in range(300):
        mid = (lo + hi) / 2
        lo, hi = (mid, hi) if f(mid) > 0 else (lo, mid)
    return lo


def inclinacao(precos):  # por bit, log2(1/eps)
    return 3 / math.log2(raiz(precos))


if __name__ == "__main__":
    print("Clifford+T: %.3f por bit (valor provado: 3)" % inclinacao([(1, 1)]))
    print("Clifford+sqrt(T), preco c da sqrt(T) em estados T:")
    for c in (1.5, 2, 3, 4, 5, 6, 8, 12):
        s = inclinacao([(1, 1), (2, c)])
        print("  c = %-4s %.3f por bit, economia %.1f%%" % (c, s, 100 * (1 - s / 3)))
    tabelas = (("(a) aplicacao direta", [1, 2.5, 3.25, 3.625, 3.8125], [10.46, 6.89, 6.05, 5.33, 4.84], [0.43, 0.22, 0.03, 0.06, 0.21]),
               ("(b) via estados magicos", [1, 3, 5, 7, 9], [10.46, 7.47, 7.19, 7.21, 7.21], [0.43, 0.15, 0.12, 0.25, 0.25]))
    for nome, C, med, er in tabelas:
        print("Mooney et al. 2021, %s (inclinacao por log10):" % nome)
        for L in range(1, 6):
            s = inclinacao([(2 ** j, C[j]) for j in range(L)]) * math.log2(10)
            if L == 1:
                base = s
            print("  Set%d: previsto %.2f (reducao %.1f%%) | medido %.2f +- %.2f (reducao %.1f%%) | medido/previsto %.3f"
                  % (L, s, 100 * (1 - s / base), med[L - 1], er[L - 1], 100 * (1 - med[L - 1] / med[0]), med[L - 1] / s))
