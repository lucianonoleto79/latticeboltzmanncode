#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
campanha_calor.py - roda os casos de Cunha & Noleto (2026) [CN26] com o lbm3d.py e monta a
comparacao com o artigo no mesmo formato (Tabelas 2 a 6 e correlacoes da Fig. 12).

  python campanha_calor.py --listar                 # casos, malhas e custo estimado
  python campanha_calor.py --rodar --processos 3    # roda o que falta (retoma do checkpoint)
  python campanha_calor.py --relatorio              # tabelas + figuras em resultados/calor_CN26/
  python acompanhar.py --todos                      # contador de todas as simulacoes em curso

Cada caso grava em casos/<tag>/ (contador em andamento.json, resultado em
resultado_calor.json).  Interromper e rodar de novo continua de onde parou.

Malhas: D (nos por diametro/altura) em tres niveis por Re, razao ~1,4 em D (~2x nos), como o
refinamento da Tabela 1 de [CN26] (~3,8x elementos por nivel em 2D com malha nao estruturada
fica entre 1,4 e 2 em D).  O passo de tempo do LBM e Dt = Dx (U fixo em unidades de rede).
"""
import argparse
import json
import os
import subprocess
import sys
import time

import numpy as np

AQUI = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, AQUI)

RES = [1, 10, 100, 200]
U_RE = {1: 0.02, 10: 0.05, 100: 0.05, 200: 0.05}           # velocidade de rede
D_NIVEIS = {1: (5, 7, 10), 10: (5, 7, 10), 100: (10, 14, 20), 200: (10, 14, 20)}  # grossa, media, fina
NIVEIS = ("Coarse", "Medium", "Fine")
FORMAS_NC = ["elipse", "placa", "hexagono", "hexagono45", "quadrado", "quadrado45"]
# metodo com solido: (rho cp)s/(rho cp)f e ks/kf.  Re = 100 e 200: propriedades ** de [CN26]
# (rho = 20 kg/m3, cp = cp,f, k = 4 W/mK).  Re = 1 e 10: [CN26] usou cobre (resfriamento de
# ~2600 D/U, inviavel aqui); o metodo nao depende das propriedades do solido desde que Bi < 0,1
# e o resfriamento seja lento perto do desenvolvimento do escoamento, entao usa-se
# Cs/Cf = 277 (Re = 1) e 104 (Re = 10), com ks/kf = 157 (Bi ~ 1e-3).
SOLIDO = {1: (277.0, 157.0), 10: (104.0, 157.0), 100: (16.9, 157.0), 200: (16.9, 157.0)}
# valores do artigo
ART_T2 = {1: (0.731, 0.775, 0.798), 10: (1.737, 1.803, 1.837), 100: (4.644, 4.987, 5.109), 200: (6.539, 7.245, 7.394)}
ART_T3 = {1: (0.726, 0.772, 0.793), 10: (1.815, 1.852, 1.863), 100: (4.330, 5.022, 5.082), 200: (5.251, 7.266, 7.414)}
ART_T5 = {1: (0.791, 0.793, 0.25, 118682, 0.798, 0.88, 8140), 10: (1.853, 1.863, 0.54, 353853, 1.837, 0.86, 4468),
          100: (5.210, 5.082, 2.46, 217085, 5.109, 1.94, 140000), 200: (7.244, 7.414, 2.35, 297635, 7.394, 2.07, 120000)}
ART_BI = {1: 2.52e-5, 10: 5.89e-5, 100: 1.66e-2, 200: 2.31e-2}
ART_T6 = {"elipse": (0.545, 1.339, 3.714, 5.210), "placa": (0.892, 1.827, 5.738, 8.516),
          "hexagono": (0.721, 1.644, 4.623, 6.493), "hexagono45": (0.793, 1.787, 4.982, 7.268),
          "quadrado": (0.631, 1.443, 3.875, 5.178), "quadrado45": (0.800, 1.770, 5.066, 7.420)}
COL_T6 = {"elipse": "Ellipse", "placa": "Vertical Plate", "hexagono": "Hexagon",
          "hexagono45": "Hexagon (tilted 45°)", "quadrado": "Square", "quadrado45": "Square (tilted 45°)"}


def casos(plano="completo"):
    """Lista de casos: dict(forma, metodo, Re, D, nivel, argumentos)."""
    lista = []
    for metodo in ("sem_solido", "com_solido"):
        for Re in RES:
            niveis = range(3) if plano == "completo" or metodo == "sem_solido" else [2]
            for n in niveis:
                lista.append(dict(forma="circulo", metodo=metodo, Re=Re, D=D_NIVEIS[Re][n], nivel=NIVEIS[n]))
    for forma in FORMAS_NC:
        for Re in RES:
            n = 2 if plano in ("completo", "reduzido") else 1
            lista.append(dict(forma=forma, metodo="sem_solido", Re=Re, D=D_NIVEIS[Re][n], nivel=NIVEIS[n]))
    for c in lista:
        a = ["--case", "calor", "--forma", c["forma"], "--metodo", c["metodo"], "--Re", "%g" % c["Re"],
             "--res", str(c["D"]), "--u-lattice", "%g" % U_RE[c["Re"]]]
        if c["metodo"] == "com_solido":
            cs, ks = SOLIDO[c["Re"]]
            a += ["--Cs-Cf", "%g" % cs, "--ks-kf", "%g" % ks]
        c["args"] = a
        c["tag"] = "calor_%s_%s_Re%g_D%d" % (c["forma"], "comsolido" if c["metodo"] == "com_solido"
                                             else "semsolido", c["Re"], c["D"])
        c["pasta"] = os.path.join(AQUI, "casos", c["tag"])
    return lista


def custo(c):
    """(nos, passos estimados) do caso."""
    import lbm3d
    Ds = lbm3d.DS_PADRAO[float(c["Re"])]
    n = (np.ceil(Ds * c["D"]) + 2) * (np.ceil(Ds * c["D"]) + 3)
    U = U_RE[c["Re"]]
    if c["metodo"] == "com_solido":
        cs = SOLIDO[c["Re"]][0]
        nu_ref = lbm3d.nu_churchill_bernstein(c["Re"])
        tdu = 3.9 * cs * c["Re"] * lbm3d.PR_AR / (4 * nu_ref) + 10
    else:
        tdu = {1: 180.0, 10: 140.0, 100: 153.0, 200: 131.0}[c["Re"]]
    return n, tdu * c["D"] / U


def resultado(c):
    p = os.path.join(c["pasta"], "resultado_calor.json")
    if os.path.exists(p):
        with open(p) as fh:
            return json.load(fh)
    return None


def listar(lista, mlups):
    tot = 0.0
    print("%-46s %-6s %10s %10s %9s  %s" % ("caso", "nivel", "nos", "passos", "horas*", "situacao"))
    for c in lista:
        n, p = custo(c)
        h = n * p / (mlups * 1e6) / 3600
        r = resultado(c)
        tot += 0 if r else h
        print("%-46s %-6s %10.0f %10.0f %9.1f  %s" % (c["tag"], c["nivel"], n, p, h,
                                                   "feito (Nu = %.4f)" % r["Nu"] if r else "-"))
    print("\n* horas de CPU a %.1f MLUPS por processo; falta ~%.0f h de CPU no total." % (mlups, tot))
    return tot


def rodar(lista, processos, max_time=0.0, menores_primeiro=False):
    py = sys.executable
    pend = [c for c in lista if resultado(c) is None]
    pend.sort(key=lambda c: (1 if menores_primeiro else -1) * np.prod(custo(c)))
    ativos = []
    print("Casos a rodar: %d (processos simultaneos: %d)" % (len(pend), processos))
    print("Contador de todas as simulacoes:  python acompanhar.py --todos")
    while pend or ativos:
        ativos = [(c, p) for c, p in ativos if p.poll() is None]
        while pend and len(ativos) < processos:
            c = pend.pop(0)
            ck = os.path.join(c["pasta"], "checkpoint", "checkpoint.npz")
            cmd = [py, os.path.join(AQUI, "lbm3d.py"), "--no-gui", "--save-every", "0",
                   "--checkpoint-every", "5000", "--report-every", "1000"]
            cmd += (["--restart", c["pasta"]] if os.path.exists(ck) else c["args"] + ["--outdir", c["pasta"]])
            if max_time:
                cmd += ["--max-time", str(max_time)]
            os.makedirs(c["pasta"], exist_ok=True)
            log = open(os.path.join(c["pasta"], "campanha.out"), "a")
            p = subprocess.Popen(cmd, cwd=AQUI, stdout=log, stderr=subprocess.STDOUT)
            print(time.strftime("%H:%M:%S"), "iniciado:", c["tag"], "(pid %d)" % p.pid, flush=True)
            ativos.append((c, p))
        time.sleep(5)
        for c, p in ativos:
            if p.poll() is not None:
                r = resultado(c)
                print(time.strftime("%H:%M:%S"), "terminou:", c["tag"],
                      "Nu = %.4f" % r["Nu"] if r else "(sem resultado: interrompido? rode de novo)", flush=True)
    print("Campanha concluida.")


# ---------------------------------------------------------------------------
# Relatorio no formato do artigo
# ---------------------------------------------------------------------------
def dif(a, b):
    return 100 * abs(a - b) / b if (a is not None and b) else float("nan")


def fmt(x, p=3):
    return "–" if x is None or (isinstance(x, float) and not np.isfinite(x)) else ("%.*f" % (p, x))


def relatorio(lista, saida):
    import lbm3d
    os.makedirs(saida, exist_ok=True)
    R = {c["tag"]: resultado(c) for c in lista}
    porchave = {(c["forma"], c["metodo"], c["Re"], c["nivel"]): R[c["tag"]] for c in lista}
    L = ["# Comparação com Cunha & Noleto (2026) [CN26]", "",
         "Resultados do lbm3d.py (LBM: cumulantes D3Q27 + energia D3Q7 TRT, 2D com nz = 1) ao lado dos "
         "valores do artigo (CVFEM). Mesma disposição das tabelas do artigo; os valores do LBM estão "
         "nas colunas marcadas **LBM**. Em vez de Coarse/Medium/Fine com malhas triangulares, o LBM usa "
         "D = nós por diâmetro (malha cartesiana uniforme).", ""]

    def tabela_malha(metodo, art, titulo):
        L.extend([titulo, "", "| Case | Re | Mesh ID | D (nós) | Nu **LBM** | dif_i,j (%) **LBM** | Nu artigo | dif_i,j (%) artigo |",
                  "|---|---|---|---|---|---|---|---|"])
        for k, Re in enumerate(RES):
            ant, anta = None, None
            for n, nv in enumerate(NIVEIS):
                r = porchave.get(("circulo", metodo, Re, nv))
                nu = r["Nu"] if r else None
                L.append("| %02d | %d | %s | %d | %s | %s | %.3f | %s |" % (
                    k + 1, Re, nv, D_NIVEIS[Re][n], fmt(nu), "–" if n == 0 or nu is None or ant is None
                    else "%.2f" % dif(ant, nu), art[Re][n], "–" if n == 0 else "%.2f" % dif(anta, art[Re][n])))
                ant, anta = nu, art[Re][n]
        L.append("")
    tabela_malha("sem_solido", ART_T2, "## Tabela 2 — Sensibilidade de malha, método sem sólido")
    tabela_malha("com_solido", ART_T3, "## Tabela 3 — Sensibilidade de malha, método com sólido")

    L.extend(["## Tabela 4 — Método com sólido (malha fina)", "",
              "| Case | Re | Nu_ref | Nu_num **LBM** | dif_num,ref **LBM** | Biot **LBM** | Nu_num artigo | dif artigo | Biot artigo | Sólido (LBM: (ρc_p)s/(ρc_p)f, ks/kf) |",
              "|---|---|---|---|---|---|---|---|---|---|"])
    for k, Re in enumerate(RES):
        r = porchave.get(("circulo", "com_solido", Re, "Fine"))
        ref = lbm3d.nu_churchill_bernstein(Re)
        a = ART_T5[Re]
        L.append("| %02d | Re = %d | %.3f | %s | %s | %s | %.3f | %.2f%% | %.2e | %g, %g |" % (
            k + 1, Re, ref, fmt(r["Nu"] if r else None), "%.2f%%" % dif(r["Nu"], ref) if r else "–",
            "%.2e" % r["Biot"] if r else "–", a[1], a[2], ART_BI[Re], *SOLIDO[Re]))
    L.append("")

    L.extend(["## Tabela 5 — Comparação dos dois métodos (malha fina)", "",
              "| Case | Re | Nu_ref | com sólido: Nu **LBM** | dif **LBM** | instantes **LBM** | sem sólido: Nu **LBM** | dif **LBM** | instantes **LBM** | com sólido: Nu artigo | dif | instantes | sem sólido: Nu artigo | dif | instantes |",
              "|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|"])
    for k, Re in enumerate(RES):
        rc = porchave.get(("circulo", "com_solido", Re, "Fine"))
        rs = porchave.get(("circulo", "sem_solido", Re, "Fine"))
        ref = lbm3d.nu_churchill_bernstein(Re)
        a = ART_T5[Re]
        L.append("| %02d | Re = %d | %.3f | %s | %s | %s | %s | %s | %s | %.3f | %.2f%% | %s | %.3f | %.2f%% | %s |" % (
            k + 1, Re, ref, fmt(rc["Nu"] if rc else None), "%.2f%%" % dif(rc["Nu"], ref) if rc else "–",
            "{:,}".format(rc["instantes"]) if rc else "–", fmt(rs["Nu"] if rs else None),
            "%.2f%%" % dif(rs["Nu"], ref) if rs else "–", "{:,}".format(rs["instantes"]) if rs else "–",
            a[1], a[2], "{:,}".format(a[3]), a[4], a[5], "{:,}".format(a[6])))
    L.append("")

    L.extend(["## Tabela 6 — Nu global para cilindros não circulares (método sem sólido)", "",
              "| Re | " + " | ".join("%s **LBM** | artigo | dif (%%)" % COL_T6[f] for f in FORMAS_NC) + " |",
              "|---|" + "---|---|---|" * len(FORMAS_NC)])
    nivel_nc = {}
    for c in lista:
        if c["forma"] != "circulo":
            nivel_nc[(c["forma"], c["Re"])] = c["nivel"]
    for k, Re in enumerate(RES):
        cel = []
        for f in FORMAS_NC:
            r = porchave.get((f, "sem_solido", Re, nivel_nc.get((f, Re), "Fine")))
            nu = r["Nu"] if r else None
            cel.append("%s | %.3f | %s" % (fmt(nu), ART_T6[f][k], "%.2f" % dif(nu, ART_T6[f][k]) if nu else "–"))
        L.append("| %d | %s |" % (Re, " | ".join(cel)))
    L.append("")

    # correlacoes Nu = C Re^m Pr^(1/3) (Fig. 12)
    L.extend(["## Fig. 12 — Correlações Nu = C·Re^m·Pr^(1/3), 1 ≤ Re ≤ 200 (mínimos quadrados)", "",
              "| Seção | C **LBM** | m **LBM** | C (Tabela 6 do artigo) | m (Tabela 6 do artigo) |",
              "|---|---|---|---|---|"])
    Pr = lbm3d.PR_AR
    corr = {}
    for f in FORMAS_NC:
        xs, ys = [], []
        for k, Re in enumerate(RES):
            r = porchave.get((f, "sem_solido", Re, nivel_nc.get((f, Re), "Fine")))
            if r:
                xs.append(Re)
                ys.append(r["Nu"])
        ca = ajuste(RES, ART_T6[f], Pr)
        cl = ajuste(xs, ys, Pr) if len(xs) >= 2 else (None, None)
        corr[f] = (cl, ca, xs, ys)
        L.append("| %s | %s | %s | %.4f | %.4f |" % (COL_T6[f], fmt(cl[0], 4), fmt(cl[1], 4), ca[0], ca[1]))
    L.extend(["", "Os coeficientes \"do artigo\" foram ajustados aqui aos valores da Tabela 6 de [CN26] com o "
              "mesmo procedimento (mínimos quadrados de Nu = C Re^m Pr^(1/3)), porque a legenda da Fig. 12 "
              "não está no texto extraído do PDF.", ""])

    # transientes: Strouhal
    L.extend(["## Escoamentos transientes (Re = 100 e 200, sem sólido, malha fina): Strouhal e C_D", "",
              "| Seção | Re | St **LBM** | C_D médio **LBM** | instantes |", "|---|---|---|---|---|"])
    for f in ["circulo"] + FORMAS_NC:
        for Re in (100, 200):
            r = porchave.get((f, "sem_solido", Re, "Fine" if f == "circulo" else nivel_nc.get((f, Re), "Fine")))
            if r:
                L.append("| %s | %d | %s | %s | %s |" % (COL_T6.get(f, "Circular"), Re, fmt(r.get("Strouhal"), 4),
                                                         fmt(r.get("C_D_medio"), 3), "{:,}".format(r["instantes"])))
    L.append("")
    with open(os.path.join(saida, "comparacao_CN26.md"), "w", encoding="utf-8") as fh:
        fh.write("\n".join(L) + "\n")
    with open(os.path.join(saida, "resultados.json"), "w") as fh:
        json.dump({c["tag"]: R[c["tag"]] for c in lista}, fh, indent=1, default=float)
    figuras_relatorio(lista, porchave, corr, nivel_nc, saida)
    print("\n".join(L))
    print("\nGravado em", saida)


def ajuste(Re, Nu, Pr):
    """C, m de Nu = C Re^m Pr^(1/3) por minimos quadrados em log."""
    x = np.log(np.asarray(Re, float))
    y = np.log(np.asarray(Nu, float) / Pr ** (1 / 3.0))
    m, lc = np.polyfit(x, y, 1)
    return float(np.exp(lc)), float(m)


def figuras_relatorio(lista, porchave, corr, nivel_nc, saida):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import lbm3d
    Pr = lbm3d.PR_AR
    # Fig. 12: correlacoes
    fig, ax = plt.subplots(figsize=(8, 5.5))
    cores = plt.rcParams["axes.prop_cycle"].by_key()["color"]
    Re = np.linspace(1, 200, 200)
    for i, f in enumerate(FORMAS_NC):
        cl, ca, xs, ys = corr[f]
        cor = cores[i % len(cores)]
        if cl[0]:
            ax.plot(Re, cl[0] * Re ** cl[1] * Pr ** (1 / 3), "-", color=cor,
                    label="%s (LBM): Nu=%.3fRe$^{%.3f}$Pr$^{1/3}$" % (COL_T6[f], cl[0], cl[1]))
            ax.plot(xs, ys, "o", color=cor)
        ax.plot(RES, ART_T6[f], "x", color=cor, ms=8)
        ax.plot(Re, ca[0] * Re ** ca[1] * Pr ** (1 / 3), ":", color=cor, lw=1)
    ax.set_xlabel("Re")
    ax.set_ylabel("Nu")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
    ax.set_title("Correlações Nu = C Re$^m$ Pr$^{1/3}$ (o, —: LBM;  x, ···: Tabela 6 de [CN26])", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "fig12_correlacoes.png"), dpi=130)
    plt.close(fig)
    # Nu x Re: LBM x artigo (circulo e nao circulares)
    fig, ax = plt.subplots(figsize=(8, 5))
    for i, f in enumerate(["circulo"] + FORMAS_NC):
        cor = cores[i % len(cores)]
        art = [ART_T5[r][4] for r in RES] if f == "circulo" else ART_T6[f]
        ax.plot(RES, art, "x--", color=cor, lw=0.8)
        nv = lambda r: "Fine" if f == "circulo" else nivel_nc.get((f, r), "Fine")
        v = [(r, porchave.get((f, "sem_solido", r, nv(r)))) for r in RES]
        v = [(r, x["Nu"]) for r, x in v if x]
        if v:
            ax.plot(*zip(*v), "o-", color=cor, label=COL_T6.get(f, "Circular"))
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("Re")
    ax.set_ylabel("Nu")
    ax.grid(alpha=0.3, which="both")
    ax.legend(fontsize=8)
    ax.set_title("Nu global: LBM (o—) x [CN26] (x--), método sem sólido", fontsize=9)
    fig.tight_layout()
    fig.savefig(os.path.join(saida, "nusselt_lbm_x_artigo.png"), dpi=130)
    plt.close(fig)
    # copia das figuras de cada caso fino (temperatura, Nu, FFT, solido)
    import shutil
    for c in lista:
        if c["nivel"] != "Fine" and not (c["forma"] != "circulo"):
            continue
        dst = os.path.join(saida, "casos", c["tag"])
        for sub, nomes in (("escoamento/final", ["temperatura.png"]),
                           ("graficos", ["nusselt_historia.png", "nusselt_7periodos.png", "nusselt_fft.png",
                                         "solido_perfil_nusselt.png", "solido_ajuste_h.png", "nusselt.csv"]),
                           ("", ["resultado_calor.json", "parametros.json", "execucao.log"])):
            for n in nomes:
                src = os.path.join(c["pasta"], sub, n)
                if os.path.exists(src):
                    os.makedirs(dst, exist_ok=True)
                    shutil.copy(src, dst)


def main():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--plano", choices=["completo", "reduzido", "rapido"], default="completo",
                   help="completo: estudo de malha nos dois metodos + nao circulares na malha fina; "
                        "reduzido: metodo com solido so na malha fina; rapido: + nao circulares na malha media")
    p.add_argument("--listar", action="store_true")
    p.add_argument("--rodar", action="store_true")
    p.add_argument("--relatorio", action="store_true")
    p.add_argument("--processos", type=int, default=3)
    p.add_argument("--mlups", type=float, default=2.5, help="velocidade por processo para a estimativa")
    p.add_argument("--max-time", type=float, default=0.0, help="limite (min) por execucao")
    p.add_argument("--menores-primeiro", action="store_true", help="roda primeiro os casos mais curtos")
    p.add_argument("--saida", default=os.path.join(AQUI, "resultados", "calor_CN26"))
    a = p.parse_args()
    lista = casos(a.plano)
    if a.listar or not (a.rodar or a.relatorio):
        listar(lista, a.mlups)
    if a.rodar:
        rodar(lista, a.processos, a.max_time, a.menores_primeiro)
    if a.relatorio:
        relatorio(lista, a.saida)


if __name__ == "__main__":
    main()
