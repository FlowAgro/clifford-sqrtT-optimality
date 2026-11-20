"""Copia tudo que o script imprime para logs/<nome>_<data_hora>.log (tela e arquivo ao mesmo tempo)."""
import os
import platform
import sys
import time

AQUI = os.path.dirname(os.path.abspath(__file__))


class _Tee:
    def __init__(self, *fluxos):
        self.fluxos = fluxos

    def write(self, s):
        for f in self.fluxos:
            f.write(s)
            f.flush()

    def flush(self):
        for f in self.fluxos:
            f.flush()


def iniciar_log(nome):
    pasta = os.path.join(AQUI, "logs")
    os.makedirs(pasta, exist_ok=True)
    caminho = os.path.join(pasta, "%s_%s.log" % (nome, time.strftime("%Y%m%d_%H%M%S")))
    arq = open(caminho, "w", encoding="utf8")
    sys.stdout = _Tee(sys.__stdout__, arq)
    sys.stderr = _Tee(sys.__stderr__, arq)
    print("# log: %s" % caminho)
    print("# comando: %s" % " ".join(sys.argv))
    print("# inicio: %s | python %s | %s" % (time.strftime("%Y-%m-%d %H:%M:%S"), platform.python_version(), platform.platform()))
    return caminho
