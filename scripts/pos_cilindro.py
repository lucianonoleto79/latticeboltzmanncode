#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""Pos-processamento do cilindro (Schaefer & Turek 3D-1Z) a partir do campo .vtk final.

No plano medio z = H/2 calcula:
  * L_a: comprimento da bolha de recirculacao, da traseira do cilindro (x = xc + R) ate o ponto da
    linha de centro (y = yc) em que u_x volta a ser >= 0 (interpolacao linear entre nos);
  * C_p(theta) = 2 (p - p_traseira) / (rho U_media^2) nos nos de fluido vizinhos a parede, com
    p = (rho - 1)/3 e p_traseira = pressao no no mais proximo de theta = 180 (assim C_p na frente
    e o Delta p frente-tras adimensional do benchmark; a pressao da saida inclui a perda de carga
    do canal); theta = 0 no ponto de estagnacao frontal;
  * angulo de separacao: angulo (a partir da frente) em que a velocidade tangencial junto a parede,
    no lado de cima, troca de sinal;
  * perfis u_x(y) em x = xc + R + {0.5, 1, 2, 4} D.
Valores fisicos usam dx = 0.1 m / D e U_media = 0.2 m/s.

Uso:
  python scripts/pos_cilindro.py casos/cylinder_D6_Re20_U0.05 [--saida <pasta>]
  (padrao: grava em <pasta do caso>/graficos/)
"""
import argparse
import glob
import json
import os
import sys

import numpy as np

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), ".."))
import lbm3d  # noqa: E402


def ler_vtk(caminho):
    """Le o VTK STRUCTURED_POINTS binario gravado por Simulation.save_vtk."""
    with open(caminho, "rb") as fh:
        dados = fh.read()
    pos = 0

    def linha():
        nonlocal pos
        fim = dados.index(b"\n", pos)
        s = dados[pos:fim].decode()
        pos = fim + 1
        return s

    while True:
        s = linha()
        if s.startswith("DIMENSIONS"):
            nx, ny, nz = map(int, s.split()[1:])
        if s.startswith("POINT_DATA"):
            break
    n = nx * ny * nz
    campos = {}
    while pos < len(dados):
        s = linha().strip()
        if not s:
            continue
        tipo, nome = s.split()[0], s.split()[1]
        if tipo == "SCALARS":
            linha()  # LOOKUP_TABLE
            a = np.frombuffer(dados, ">f4", n, pos)
            pos += 4 * n
            campos[nome] = a.reshape(nz, ny, nx).transpose(2, 1, 0).astype(float)
        elif tipo == "VECTORS":
            a = np.frombuffer(dados, ">f4", 3 * n, pos)
            pos += 12 * n
            campos[nome] = a.reshape(nz, ny, nx, 3).transpose(3, 2, 1, 0).astype(float)
    return campos


def montar_caso(pasta):
    par = json.load(open(os.path.join(pasta, "parametros.json")))
    a = lbm3d.parse_args([])
    for k, v in par["argumentos"].items():
        setattr(a, k, v)
    return lbm3d.Cylinder(a)


def analisar(pasta):
    vtks = sorted(glob.glob(os.path.join(pasta, "escoamento", "final", "*.vtk")))
    if not vtks:
        raise SystemExit("sem .vtk em %s/escoamento/final (o .vtk nao e versionado: use a pasta em casos/)"
                         % pasta)
    case = montar_caso(pasta)
    c = ler_vtk(vtks[-1])
    rho, u, solido = c["density"], c["velocity"], c["solid"] > 0.5
    k = case.H // 2 + 1                      # plano z ~ H/2 (no k fica em z = k - 0.5)
    ux, uy = u[0][:, :, k], u[1][:, :, k]
    p = (rho[:, :, k] - 1.0) / 3.0
    sol = solido[:, :, k]
    dx, R, xc, yc, D = case.dx_fisico, case.R, case.xc, case.yc, case.D
    Ub = case.Ubar

    # --- L_a na linha y = yc (interpolacao entre as duas linhas de nos vizinhas)
    j0 = int(np.floor(yc + 0.5))
    w = yc + 0.5 - j0
    linha_ux = (1 - w) * ux[:, j0] + w * ux[:, j0 + 1]
    xs = np.arange(ux.shape[0]) - 0.5
    atras = xs > xc + R
    xa, va = xs[atras], linha_ux[atras]
    L_a = np.nan
    neg = np.where(va < 0)[0]
    if neg.size:
        i = neg[-1]
        if i + 1 < va.size:
            x0 = xa[i] - va[i] * (xa[i + 1] - xa[i]) / (va[i + 1] - va[i])
            L_a = x0 - (xc + R)

    # --- C_p e velocidade tangencial nos nos de fluido vizinhos ao cilindro
    ii, jj = np.nonzero(~sol)
    X, Y = ii - 0.5 - xc, jj - 0.5 - yc
    r = np.hypot(X, Y)
    perto = (r > R) & (r < R + 1.5) & (ii > 1) & (ii < sol.shape[0] - 2)
    th = np.degrees(np.arctan2(Y[perto], -X[perto]))  # 0 = frente, +90 = topo (y > yc)
    cp = 2.0 * p[ii[perto], jj[perto]] / (Ub ** 2)
    # tangencial no sentido frente -> tras pelo lado de cima: t = (sin(phi), -cos(phi))... com phi
    # medido a partir de -x: x = -r cos, y = r sin  ->  t = d/dtheta = (r sin, r cos)/r
    ang = np.radians(th)
    ut = ux[ii[perto], jj[perto]] * np.sin(ang) + uy[ii[perto], jj[perto]] * np.cos(ang)
    ordem = np.argsort(th)
    th, cp, ut = th[ordem], cp[ordem], ut[ordem]
    cp_saida_tras = cp[np.argmin(np.abs(np.abs(th) - 180))]
    cp = cp - cp_saida_tras  # referencia: traseira
    topo = th > 0
    sep = np.nan
    t_top, u_top = th[topo], ut[topo]
    troca = np.where((u_top[:-1] > 0) & (u_top[1:] <= 0))[0]
    if troca.size:
        i = troca[-1]
        sep = t_top[i] + (t_top[i + 1] - t_top[i]) * u_top[i] / (u_top[i] - u_top[i + 1])

    # --- perfis
    perfis = {}
    for f in (0.5, 1.0, 2.0, 4.0):
        x = xc + R + f * D
        i = int(np.floor(x + 0.5))
        wx = x + 0.5 - i
        perfis[f] = (1 - wx) * ux[i, :] + wx * ux[i + 1, :]

    return dict(case=case, L_a=L_a, sep=sep, th=th, cp=cp, ut=ut, perfis=perfis, k=k,
                linha=(xa, va), vtk=vtks[-1], cp_frente=cp[np.argmin(np.abs(th))],
                cp_tras_saida=cp_saida_tras)


def main():
    ap = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    ap.add_argument("pasta")
    ap.add_argument("--saida", default="")
    a = ap.parse_args()
    r = analisar(a.pasta)
    case = r["case"]
    saida = a.saida or os.path.join(a.pasta, "graficos")
    os.makedirs(saida, exist_ok=True)
    dx, D = case.dx_fisico, case.D
    txt = ["Pos-processamento do cilindro - %s (campo %s, plano z = no %d)"
           % (a.pasta, os.path.basename(r["vtk"]), r["k"]),
           "  L_a (bolha de recirculacao) = %.3f nos = %.4f m = %.3f D"
           % (r["L_a"], r["L_a"] * dx, r["L_a"] / D),
           "  angulo de separacao (lado de cima, a partir da frente) = %.1f graus" % r["sep"],
           "  C_p no ponto de estagnacao frontal (no mais proximo) = %.3f" % r["cp_frente"],
           "  (C_p = 2 (p - p_traseira) / (rho U_media^2), nos de fluido a menos de 1.5 dx da parede;",
           "   em relacao a saida, a traseira tem C_p = %.3f)" % r["cp_tras_saida"],
           "  Com D = %d ha %d nos no raio: L_a e o angulo de separacao sao estimativas grosseiras."
           % (D, int(case.R))]
    txt = "\n".join(txt) + "\n"
    open(os.path.join(saida, "pos_cilindro.txt"), "w", encoding="utf-8").write(txt)
    print(txt)

    from matplotlib.figure import Figure
    fig = Figure(figsize=(14, 4.4))
    ax = fig.add_subplot(1, 3, 1)
    ax.plot(r["th"], r["cp"], "o", ms=3)
    ax.set_xlabel("theta [graus] (0 = frente, +90 = topo)")
    ax.set_ylabel("C_p (ref.: traseira)")
    ax.set_title("Pressao junto a parede")
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(1, 3, 2)
    xa, va = r["linha"]
    ax.plot((xa - case.xc - case.R) / D, va / case.Ubar, "k.-", ms=3)
    ax.axhline(0, color="0.5", lw=0.8)
    if np.isfinite(r["L_a"]):
        ax.axvline(r["L_a"] / D, color="r", ls="--", label="L_a = %.2f D" % (r["L_a"] / D))
        ax.legend(fontsize=8)
    ax.set_xlim(0, 6)
    ax.set_xlabel("(x - x_traseira) / D")
    ax.set_ylabel("u_x / U_media na linha de centro")
    ax.set_title("Bolha de recirculacao")
    ax.grid(alpha=0.3)
    ax = fig.add_subplot(1, 3, 3)
    y = (np.arange(case.shape[1]) - 0.5) * dx
    for f, perfil in r["perfis"].items():
        ax.plot(perfil / case.Ubar, y, label="x = traseira + %g D" % f)
    ax.set_xlabel("u_x / U_media")
    ax.set_ylabel("y [m]")
    ax.set_title("Perfis de velocidade (z = H/2)")
    ax.legend(fontsize=7)
    ax.grid(alpha=0.3)
    fig.suptitle("Cilindro 3D-1Z - %s" % os.path.basename(os.path.normpath(a.pasta)))
    fig.tight_layout()
    out = os.path.join(saida, "pos_cilindro.png")
    fig.savefig(out, dpi=120)
    print("grafico:", out)


if __name__ == "__main__":
    main()
