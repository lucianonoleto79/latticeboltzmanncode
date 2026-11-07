#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Resumo de uma execucao do lbm3d.py para a pagina da execucao no GitHub Actions.

Uso: python scripts/resumo_execucao.py <pasta_saida> [ID_DA_EXECUCAO]

Imprime Markdown (o fluxo acrescenta em $GITHUB_STEP_SUMMARY) e, se a variavel
GITHUB_OUTPUT existir, grava terminou=true|false para o fluxo decidir se continua.
"""
import json
import os
import sys


def main():
    pasta = sys.argv[1]
    run_id = sys.argv[2] if len(sys.argv) > 2 else ""
    caminho = os.path.join(pasta, "estado.json")
    if not os.path.exists(caminho):
        print("## Simulacao LBM\n\nA simulacao nao gravou `estado.json` (erro antes de comecar?). "
              "Veja o log da etapa *Rodar simulacao*.")
        terminou = True
    else:
        e = json.load(open(caminho))
        terminou = e["terminou"] or e["divergiu"]
        if e["divergiu"]:
            situacao = "**divergiu** (NaN)"
        elif e["convergiu"]:
            situacao = "**convergiu** (residuo %.3e < tol %.1e)" % (e["residuo"], e["tol"])
        elif e["terminou"]:
            situacao = "terminou (max-iter atingido; residuo %.3e)" % e["residuo"]
        else:
            situacao = "**interrompida** pelo limite de tempo (residuo %.3e)" % e["residuo"]
        print("## Simulacao LBM: caso `%s`\n" % e["caso"])
        print("| iteracao | situacao |\n|---|---|\n| %d | %s |\n" % (e["it"], situacao))
        val = os.path.join(pasta, "validacao.txt")
        if os.path.exists(val):
            print("### Validacao\n\n```\n%s```\n" % open(val).read())
        if not terminou:
            print("O checkpoint esta no artefato `resultado`. Para continuar, rode o fluxo "
                  "*Simulacao LBM* com **reiniciar_de = `%s`** (ou deixe *continuar* marcado "
                  "para isso acontecer sozinho)." % run_id)
    print("\nBaixe o artefato `resultado` desta execucao e abra `convergencia.html` no "
          "navegador para ver a convergencia; o `.vtk` abre no ParaView.")
    if os.environ.get("GITHUB_OUTPUT"):
        with open(os.environ["GITHUB_OUTPUT"], "a") as fh:
            fh.write("terminou=%s\n" % str(terminou).lower())


if __name__ == "__main__":
    main()
