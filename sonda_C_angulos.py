"""SONDA C: em circuitos publicos, quantos angulos de rotacao DISTINTOS existem por circuito?

IDEIA. Se poucos angulos se repetem muito, da para pagar a sintese lenta (60 s ou certificada) uma vez por
angulo e reaproveitar: a folga de 1-3% do cyclosynth vira gratis em tempo de compilacao.

DADOS. QASMBench (github.com/pnnl/QASMBench), clonado em ~/dados_externos/QASMBench (pasta propria; os .qasm
sao lidos so como texto, nada deles e executado). Angulos sao avaliados por um analisador restrito
(numeros, pi, + - * / e parenteses).

DEFINICOES.
  Rotacao que precisa de sintese = rz/rx/ry/p/u1/crz/cry/crx/cp/cu1 com angulo que NAO e multiplo de pi/4,
  e u2/u3/u/cu3 com algum angulo assim (a tripla conta como um alvo).
  R = numero dessas rotacoes; D = numero de alvos distintos (angulo modulo 2*pi, arredondado a 1e-9).

REGRA FIXADA ANTES DE RODAR (19/11/2026):
  C e viavel se, entre os circuitos com R >= 1000, a mediana de D/R for <= 0,10.
  Descritivo: o mesmo para os nomes que lembram simulacao (vqe, qaoa, ising, ham, trotter, heis, h2, lih, qpe, hhl...).

USO (WSL ou qualquer terminal com git e python):
  python sonda_C_angulos.py               # clona (se preciso) e mede
  python sonda_C_angulos.py --pasta X     # usa outra pasta de .qasm
Log automatico em logs/sonda_C_<data>.log e tabela em logs/sonda_C_<data>.csv.
"""
import argparse
import csv
import math
import os
import re
import statistics
import subprocess

import registro

REPO = "https://github.com/pnnl/QASMBench.git"
UM = re.compile(r"^\s*(rz|rx|ry|p|u1|crz|cry|crx|cp|cu1)\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)", re.I)
TRES = re.compile(r"^\s*(u2|u3|u|cu3)\s*\(([^()]*(?:\([^()]*\)[^()]*)*)\)", re.I)
PERMITIDO = re.compile(r"^[0-9eE.+\-*/() pi]*$")
SIM = re.compile(r"vqe|qaoa|ising|ham|trotter|heis|h2|lih|qpe|hhl|chem|mol|hubbard|qft", re.I)


def avaliar(txt):
    txt = txt.strip()
    if not PERMITIDO.match(txt) or "**" in txt:
        raise ValueError(txt)
    return float(eval(txt.replace("pi", "(%r)" % math.pi), {"__builtins__": {}}, {}))


def trivial(x):
    k = x / (math.pi / 4)
    return abs(k - round(k)) < 1e-9


def chave(x):
    return round(math.fmod(x, 2 * math.pi) % (2 * math.pi), 9)


def medir(arq):
    R, alvos, ruins = 0, set(), 0
    for linha in open(arq, encoding="utf8", errors="replace"):
        m = UM.match(linha) or TRES.match(linha)
        if not m:
            continue
        try:
            angs = [avaliar(s) for s in m.group(2).split(",")]
        except Exception:
            ruins += 1
            continue
        if all(trivial(a) for a in angs):
            continue
        R += 1
        alvos.add((m.group(1).lower() in ("rx", "crx"), m.group(1).lower() in ("ry", "cry")) + tuple(chave(a) for a in angs))
    return R, len(alvos), ruins


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pasta", default=os.path.expanduser("~/dados_externos/QASMBench"))
    args = p.parse_args()
    log = registro.iniciar_log("sonda_C")
    if not os.path.isdir(args.pasta):
        os.makedirs(os.path.dirname(args.pasta), exist_ok=True)
        print("clonando %s em %s ..." % (REPO, args.pasta), flush=True)
        subprocess.run(["git", "clone", "--depth", "1", REPO, args.pasta], check=True)
    rev = subprocess.run(["git", "-C", args.pasta, "rev-parse", "--short", "HEAD"], capture_output=True, text=True).stdout.strip()
    print("dados: %s (commit %s)" % (args.pasta, rev or "?"))

    linhas = []
    for raiz, _, nomes in os.walk(args.pasta):
        for nm in nomes:
            if nm.endswith(".qasm"):
                arq = os.path.join(raiz, nm)
                R, D, ruins = medir(arq)
                if R:
                    linhas.append((os.path.relpath(arq, args.pasta), R, D, ruins))
    csv_out = log[:-4] + ".csv"
    with open(csv_out, "w", encoding="utf8", newline="") as f:
        w = csv.writer(f)
        w.writerow(["arquivo", "R_rotacoes", "D_distintos", "D_sobre_R", "linhas_nao_lidas"])
        for a, R, D, ru in sorted(linhas, key=lambda x: -x[1]):
            w.writerow([a, R, D, "%.4f" % (D / R), ru])

    print("circuitos com rotacoes nao triviais: %d | linhas de angulo nao lidas: %d" % (len(linhas), sum(x[3] for x in linhas)))
    print("\n20 maiores (R = rotacoes a sintetizar, D = alvos distintos):")
    for a, R, D, _ in sorted(linhas, key=lambda x: -x[1])[:20]:
        print("  %-60s R=%9d D=%8d D/R=%.4f" % (a[:60], R, D, D / R))

    def resumo(nome, sel):
        if not sel:
            print("%s: nenhum circuito" % nome)
            return None
        razoes = [D / R for _, R, D, _ in sel]
        tot_R, tot_D = sum(x[1] for x in sel), sum(x[2] for x in sel)
        med = statistics.median(razoes)
        print("%s: n = %d | mediana D/R = %.4f | soma R = %d, soma D = %d (D/R agregado %.4f)"
              % (nome, len(sel), med, tot_R, tot_D, tot_D / tot_R))
        return med

    print()
    med = resumo("R >= 1000", [x for x in linhas if x[1] >= 1000])
    resumo("R >= 1000 e nome de simulacao", [x for x in linhas if x[1] >= 1000 and SIM.search(x[0])])
    resumo("todos", linhas)
    print("\nVEREDITO C:", "sem dados" if med is None else
          ("viavel: poucos angulos distintos, sintese lenta amortizada" if med <= 0.10 else
           "nao viavel: a maioria dos angulos e distinta"))
    print("tabela por circuito em:", csv_out)


if __name__ == "__main__":
    main()
