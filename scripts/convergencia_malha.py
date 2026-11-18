#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Analise de convergencia de malha do cilindro (Schaefer & Turek 3D-1Z).

Le as pastas das rodadas (padrao: casos/ e resultados/ com cylinder_D*_Re20_U0.05), pega os
valores finais de C_D, C_L e dp (ultima linha de graficos/convergencia.csv) e grava em --saida:
  convergencia_malha.png  grandezas x 1/D com as faixas de referencia e a extrapolacao
  convergencia_malha.txt  tabela, ordem observada e extrapolacao de Richardson

Extrapolacao (tres malhas mais finas, razao de refino nao constante):
  f(h) = f0 + C h^p,  h = 1/D;  p resolve (f1 - f2)/(f2 - f3) = (h1^p - h2^p)/(h2^p - h3^p)
  (so quando as diferencas tem o mesmo sinal, isto e, convergencia monotonica).

Uso:
  python scripts/convergencia_malha.py
  python scripts/convergencia_malha.py --sufixo U0.05 --saida resultados/cylinder_convergencia_malha
"""
import argparse
import glob
import json
import os
import re
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import acompanhar  # noqa: E402
from lbm3d import Cylinder  # noqa: E402

NOMES = {"C_D": "C_D", "C_L": "C_L", "dp": "Δp [Pa]"}


def coletar(sufixo):
    rodadas = {}
    for base in ("resultados", "casos"):  # casos/ (rodada local) tem prioridade
        for pasta in glob.glob(os.path.join(base, "cylinder_D*_Re20_%s" % sufixo)):
            m = re.search(r"cylinder_D(\d+)_", pasta)
            meta, c = acompanhar.ler(pasta)
            if not m or c is None:
                continue
            est = {}
            if os.path.exists(os.path.join(pasta, "estado.json")):
                est = json.load(open(os.path.join(pasta, "estado.json")))
            rodadas[int(m.group(1))] = dict(
                pasta=pasta, it=int(c["it"][-1]), res=float(c["res"][-1]),
                convergiu=bool(est.get("convergiu", False)),
                **{k: float(c[k][-1]) for k in NOMES})
    return dict(sorted(rodadas.items()))


def richardson(h, f):
    """h e f das tres malhas mais finas, h decrescente.  Devolve (p, f0) ou None."""
    (h1, h2, h3), (f1, f2, f3) = h, f
    d12, d23 = f1 - f2, f2 - f3
    if d12 == 0 or d23 == 0 or np.sign(d12) != np.sign(d23):
        return None
    alvo = d12 / d23
    g = lambda p: (h1 ** p - h2 ** p) / (h2 ** p - h3 ** p) - alvo
    lo, hi = 0.05, 8.0
    if g(lo) * g(hi) > 0:
        return None
    for _ in range(200):
        mid = 0.5 * (lo + hi)
        if g(lo) * g(mid) <= 0:
            hi = mid
        else:
            lo = mid
    p = 0.5 * (lo + hi)
    C = d23 / (h2 ** p - h3 ** p)
    return p, f3 - C * h3 ** p


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("--sufixo", default="U0.05", help="sufixo das pastas (velocidade de rede)")
    ap.add_argument("--saida", default="resultados/cylinder_convergencia_malha")
    a = ap.parse_args()
    rod = coletar(a.sufixo)
    if not rod:
        sys.exit("nenhuma rodada cylinder_D*_Re20_%s encontrada em casos/ ou resultados/" % a.sufixo)
    os.makedirs(a.saida, exist_ok=True)
    Ds = np.array(list(rod))
    h = 1.0 / Ds
    linhas = ["Convergencia de malha - cilindro Schaefer & Turek 3D-1Z, Re = 20, K17, Um = %s"
              % a.sufixo.lstrip("U"), "",
              "%4s %8s %10s %10s %10s %10s  %s" % ("D", "it", "residuo", "C_D", "C_L", "dp [Pa]", "pasta")]
    for D, r in rod.items():
        linhas.append("%4d %8d %10.2e %10.5f %10.5f %10.5f  %s%s" % (
            D, r["it"], r["res"], r["C_D"], r["C_L"], r["dp"], r["pasta"],
            "" if r["convergiu"] else "  (NAO CONVERGIU: valor provisorio)"))
    linhas += ["", "Faixas de referencia: " + ", ".join(
        "%s [%g, %g]" % (k, *Cylinder.REF[k]) for k in NOMES)]
    extrap = {}
    if len(Ds) >= 3:
        linhas += ["", "Extrapolacao de Richardson (malhas D = %s):" % ", ".join(map(str, Ds[-3:]))]
        for k in NOMES:
            f = [rod[D][k] for D in Ds[-3:]]
            r = richardson(h[-3:], f)
            if r is None:
                linhas.append("  %-4s convergencia nao monotonica: sem extrapolacao" % k)
            else:
                p, f0 = r
                extrap[k] = f0
                lo, hi = Cylinder.REF[k]
                linhas.append("  %-4s ordem observada p = %.2f   valor extrapolado = %.5f  (%s)"
                              % (k, p, f0, "dentro da faixa" if lo <= f0 <= hi else "fora da faixa"))
    else:
        linhas += ["", "Extrapolacao de Richardson: precisa de pelo menos 3 malhas."]
    txt = "\n".join(linhas) + "\n"
    open(os.path.join(a.saida, "convergencia_malha.txt"), "w", encoding="utf-8").write(txt)
    print(txt)

    from matplotlib.figure import Figure
    fig = Figure(figsize=(13, 4.2))
    for n, (k, nome) in enumerate(NOMES.items()):
        ax = fig.add_subplot(1, 3, n + 1)
        lo, hi = Cylinder.REF[k]
        ax.axhspan(lo, hi, color="#1baf7a", alpha=0.18, label="faixa de referencia")
        v = [rod[D][k] for D in Ds]
        ax.plot(h, v, "ko-", label="LBM K17")
        for D, x, y in zip(Ds, h, v):
            ax.annotate("D=%d" % D, (x, y), textcoords="offset points", xytext=(4, 4), fontsize=7)
        if k in extrap:
            ax.plot([0], [extrap[k]], "r*", ms=10, label="Richardson (D→∞)")
        ax.set_xlim(left=0)
        ax.set_xlabel("1/D")
        ax.set_title(nome)
        ax.grid(alpha=0.3)
        ax.legend(fontsize=7)
    fig.suptitle("Convergencia de malha - cilindro 3D-1Z, Re = 20, Um = %s" % a.sufixo.lstrip("U"))
    fig.tight_layout()
    out = os.path.join(a.saida, "convergencia_malha.png")
    fig.savefig(out, dpi=120)
    print("grafico:", out)


if __name__ == "__main__":
    main()
