#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Compara rodadas do cilindro (Schaefer & Turek 3D-1Z) entre si: modelos de colisao, velocidade
de rede (Mach), limitador etc.  Le o valor final de C_D, C_L e dp (ultima linha de
graficos/convergencia.csv) de cada pasta e grava em --saida:
  <nome>.png  barras de C_D, C_L e dp com as faixas de referencia
  <nome>.txt  tabela com o desvio em relacao a faixa

Uso:
  python scripts/comparar_cilindro.py --nome comparacao_modelos \\
      --saida resultados/cylinder_comparacao_modelos \\
      "K17=casos/cylinder_D6_Re20_U0.05" "cumulante=casos/cylinder_D6_Re20_U0.05_cumulant" ...
"""
import argparse
import json
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import acompanhar  # noqa: E402
from lbm3d import Cylinder  # noqa: E402

NOMES = {"C_D": "C_D", "C_L": "C_L", "dp": "Δp [Pa]"}


def desvio(v, lo, hi):
    """Distancia relativa ate a faixa (0 se dentro), em % do limite mais proximo."""
    if lo <= v <= hi:
        return 0.0
    ref = lo if v < lo else hi
    return 100.0 * (v - ref) / ref


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("rodadas", nargs="+", help="ROTULO=PASTA")
    ap.add_argument("--nome", required=True)
    ap.add_argument("--saida", required=True)
    ap.add_argument("--titulo", default="")
    a = ap.parse_args()
    dados = []
    for item in a.rodadas:
        rot, pasta = item.rsplit("=", 1)
        meta, c = acompanhar.ler(pasta)
        if c is None:
            sys.exit("sem convergencia.csv em %s" % pasta)
        est = json.load(open(os.path.join(pasta, "estado.json")))
        dados.append(dict(rot=rot, pasta=pasta, it=int(c["it"][-1]), res=float(c["res"][-1]),
                          conv=est["convergiu"], div=est["divergiu"],
                          **{k: float(c[k][-1]) for k in NOMES}))
    os.makedirs(a.saida, exist_ok=True)
    L = [a.titulo or a.nome, "",
         "%-14s %7s %9s %10s %9s %10s %9s %10s %9s  %s" % (
             "rodada", "it", "residuo", "C_D", "desvio%", "C_L", "desvio%", "dp [Pa]", "desvio%", "pasta")]
    for d in dados:
        dv = {k: desvio(d[k], *Cylinder.REF[k]) for k in NOMES}
        L.append("%-14s %7d %9.2e %10.5f %+9.2f %10.5f %+9.1f %10.5f %+9.2f  %s%s" % (
            d["rot"], d["it"], d["res"], d["C_D"], dv["C_D"], d["C_L"], dv["C_L"], d["dp"], dv["dp"],
            d["pasta"], "" if d["conv"] else ("  (DIVERGIU)" if d["div"] else "  (NAO CONVERGIU)")))
    L += ["", "desvio% = distancia ate a faixa de referencia, em % do limite mais proximo (0 = dentro).",
          "Faixas: " + ", ".join("%s [%g, %g]" % (k, *Cylinder.REF[k]) for k in NOMES)]
    txt = "\n".join(L) + "\n"
    open(os.path.join(a.saida, a.nome + ".txt"), "w", encoding="utf-8").write(txt)
    print(txt)

    from matplotlib.figure import Figure
    fig = Figure(figsize=(13, 4.2))
    x = range(len(dados))
    for n, (k, nome) in enumerate(NOMES.items()):
        ax = fig.add_subplot(1, 3, n + 1)
        lo, hi = Cylinder.REF[k]
        ax.axhspan(lo, hi, color="#1baf7a", alpha=0.25, label="faixa de referencia")
        v = [d[k] for d in dados]
        ax.bar(x, v, color="#2a78d6", width=0.6)
        for i, y in zip(x, v):
            ax.annotate("%.4g" % y, (i, y), ha="center", va="bottom", fontsize=8)
        ax.set_xticks(list(x))
        ax.set_xticklabels([d["rot"] for d in dados], fontsize=8)
        top = max(max(v), hi)
        bot = min(min(v), lo)
        pad = 0.15 * (top - bot) if top > bot else 0.1 * abs(top)
        ax.set_ylim(bot - pad if k != "C_L" else 0, top + pad)
        ax.set_title(nome)
        ax.grid(alpha=0.3, axis="y")
        ax.legend(fontsize=7, loc="lower right")
    fig.suptitle(a.titulo or a.nome)
    fig.tight_layout()
    out = os.path.join(a.saida, a.nome + ".png")
    fig.savefig(out, dpi=120)
    print("grafico:", out)


if __name__ == "__main__":
    main()
