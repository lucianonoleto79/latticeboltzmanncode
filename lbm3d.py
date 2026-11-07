#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
lbm3d.py - Escoamentos laminares 3D com o metodo de lattice-Boltzmann (LBM)
============================================================================

Metodo numerico: LBM D3Q27 com operador de colisao CUMULANTE PARAMETRIZADO
("K17", difusao com precisao de 4a ordem):

  [G15] M. Geier, M. Schoenherr, A. Pasquali, M. Krafczyk, "The cumulant lattice
        Boltzmann equation in three dimensions: Theory and validation",
        Computers & Mathematics with Applications 70 (2015) 507-547.
  [G17] M. Geier, A. Pasquali, M. Schoenherr, "Parametrization of the cumulant
        lattice Boltzmann method for fourth order accurate diffusion. Part I:
        Derivation and validation", J. Comput. Phys. 348 (2017) 862-888.

Condicoes de contorno:
  [L94] A.J.C. Ladd, J. Fluid Mech. 271 (1994) 285-309  (bounce-back, parede movel)
  [B01] M. Bouzidi, M. Firdaouss, P. Lallemand, Phys. Fluids 13 (2001) 3452
        (bounce-back interpolado linear para paredes curvas)
  [G02] Z. Guo, C. Zheng, B. Shi, Chinese Physics 11 (2002) 366
        (extrapolacao de nao-equilibrio, saida com pressao imposta)
  [M02] R. Mei, D. Yu, W. Shyy, L.-S. Luo, Phys. Rev. E 65 (2002) 041203
        (forca por troca de quantidade de movimento)

Casos de teste / validacao:
  shearwave : onda de cisalhamento com decaimento viscoso (solucao exata) [G17]
  duct      : escoamento em duto retangular (solucao em serie, F.M. White,
              Viscous Fluid Flow, 3a ed., 2006, Eq. 3-48)
  cylinder  : benchmark 3D-1Z (Re=20) de M. Schaefer, S. Turek, "Benchmark
              computations of laminar flow around a cylinder", Notes on
              Numerical Fluid Mechanics 52 (1996) 547-566
  spheres   : escoamento de Stokes num arranjo cubico simples de esferas, arraste
              comparado com A.S. Sangani, A. Acrivos, Int. J. Multiphase Flow 8 (1982)
              343-360 (valores tabelados em Bogner et al., IJMF 68 (2015) 71-79)
  cavity    : cavidade cubica com tampa movel (demonstracao)

Acompanhamento da convergencia: a cada 100 it. uma linha e acrescentada a
<saida>/convergencia.csv e o painel <saida>/convergencia.html e regravado:
abra-o no navegador (recarrega sozinho) ou use  python acompanhar.py <saida> --servir

Controles (janela grafica):
  botao/tecla  P  : pausar / continuar
  botao/tecla  S  : salvar resultados parciais agora (VTK + checkpoint)
  caixa "Auto-salvar": salva automaticamente a cada 100 iteracoes
  botao/tecla  Q  : salvar checkpoint e sair (reinicie com --restart)
Sem janela (--no-gui): Ctrl+C salva o checkpoint e sai.

Exemplos:
  python lbm3d.py --case duct
  python lbm3d.py --case cylinder --res 10
  python lbm3d.py --restart saida_cylinder/checkpoint.npz
  python lbm3d.py --case spheres --res 48 --no-gui --max-time 320   # para em 320 min
"""

import argparse
import json
import os
import signal
import sys
import time

import numpy as np

# ---------------------------------------------------------------------------
# Rede D3Q27.  As populacoes sao guardadas como f[a, b, c, x, y, z], com a
# velocidade discreta (a-1, b-1, c-1).  Pesos = produto dos pesos 1D
# (1/6, 2/3, 1/6), cs^2 = 1/3.
# ---------------------------------------------------------------------------
W1 = np.array([1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0])
W = W1[:, None, None] * W1[None, :, None] * W1[None, None, :]
CV = np.array([(a - 1, b - 1, c - 1) for a in range(3) for b in range(3) for c in range(3)])
OPP = 26 - np.arange(27)  # indice da direcao oposta (d = 9a + 3b + c)

# Momentos centrais (normalizados por rho) da Maxwelliana discreta:
# k_abc^eq = e_a e_b e_c, e = (1, 0, 1/3)
E1 = np.array([1.0, 0.0, 1.0 / 3.0])
KEQ = E1[:, None, None] * E1[None, :, None] * E1[None, None, :]


# ---------------------------------------------------------------------------
# Transformacoes "chimera" para momentos centrais [G15, Eqs. (J.1)-(J.6)]
#   k0 = f- + f0 + f+
#   k1 = (f+ - f-) - u k0
#   k2 = (f+ + f-) - 2u (f+ - f-) + u^2 k0
# e a inversa
#   f- = ((u^2 - u) k0 + (2u - 1) k1 + k2) / 2
#   f0 = (1 - u^2) k0 - 2u k1 - k2
#   f+ = ((u^2 + u) k0 + (2u + 1) k1 + k2) / 2
# ---------------------------------------------------------------------------
def _fwd(m, u, axis):
    v = np.moveaxis(m, axis, 0)
    fm, f0, fp = v[0], v[1], v[2]
    k0 = fm + f0 + fp
    s = fp - fm
    k2 = (fp + fm) - 2.0 * u * s + u * u * k0
    k1 = s - u * k0
    v[0] = k0
    v[1] = k1
    v[2] = k2


def _bwd(m, u, axis):
    v = np.moveaxis(m, axis, 0)
    k0, k1, k2 = v[0], v[1], v[2]
    u2 = u * u
    fm = 0.5 * ((u2 - u) * k0 + (2.0 * u - 1.0) * k1 + k2)
    fp = 0.5 * ((u2 + u) * k0 + (2.0 * u + 1.0) * k1 + k2)
    f0 = (1.0 - u2) * k0 - 2.0 * u * k1 - k2
    v[0] = fm
    v[1] = f0
    v[2] = fp


def equilibrium(rho, ux, uy, uz):
    """Equilibrio D3Q27 (momentos centrais da Maxwelliana)."""
    m = np.empty((3, 3, 3) + rho.shape)
    m[...] = KEQ.reshape((3, 3, 3) + (1,) * rho.ndim)
    _bwd(m, ux, 0)
    _bwd(m, uy, 1)
    _bwd(m, uz, 2)
    m *= rho
    return m


def macroscopic(f, force=(0.0, 0.0, 0.0)):
    rho = f.sum(axis=(0, 1, 2))
    jx = f[2].sum(axis=(0, 1)) - f[0].sum(axis=(0, 1))
    jy = f[:, 2].sum(axis=(0, 1)) - f[:, 0].sum(axis=(0, 1))
    jz = f[:, :, 2].sum(axis=(0, 1)) - f[:, :, 0].sum(axis=(0, 1))
    return (rho, (jx + 0.5 * force[0]) / rho, (jy + 0.5 * force[1]) / rho,
            (jz + 0.5 * force[2]) / rho)


# ---------------------------------------------------------------------------
# Taxas de relaxacao do cumulante parametrizado [G17], com omega2 = 1
# (formulas de [G17, Eqs. (14)-(18)] avaliadas em omega2 = 1)
# ---------------------------------------------------------------------------
def k17_parameters(w):
    w3 = 8.0 * (w - 2.0) * (1.0 + 2.0 * w) / (7.0 * w * w - 14.0 * w - 8.0)
    w4 = 8.0 * (w - 2.0) * (4.0 * w - 7.0) / (9.0 * w * w - 50.0 * w + 56.0)
    w5 = (24.0 * (w - 2.0) * (3.0 * w * w - 7.0 * w - 2.0)
          / (29.0 * w ** 3 - 130.0 * w * w + 152.0 * w + 48.0))
    A = (4.0 + 2.0 * w - 3.0 * w * w) / (2.0 - 7.0 * w + 5.0 * w * w)
    B = (4.0 + 28.0 * w - 14.0 * w * w) / (6.0 - 21.0 * w + 15.0 * w * w)
    return dict(w2=1.0, w3=w3, w4=w4, w5=w5, A=A, B=B)


OMEGA_MIN_K17 = 1.75


def collision_parameters(omega, model, limiter=None):
    """Monta o dicionario de parametros da colisao."""
    prm = dict(model=model, w1=omega, limiter=limiter)
    if model == "k17":
        # Com omega2 = 1, A tem polo em omega = 1 e w4 tem polo em omega ~ 1.556
        # (w4 fica fora de (0,2) para 1.43 < omega < 1.75).  Nos testes de onda de
        # cisalhamento o K17 so e estavel e mais preciso para omega >= 1.75
        # (nu <= 0.024), que e o regime de baixa viscosidade para o qual foi criado.
        # Abaixo disso usa-se o cumulante classico [G15].
        if omega < OMEGA_MIN_K17:
            print("AVISO: omega = %.3f < %.2f (nu = %.4f > %.4f): fora da faixa da "
                  "parametrizacao K17 [G17]; usando cumulante classico [G15]."
                  % (omega, OMEGA_MIN_K17, (1 / omega - 0.5) / 3, (1 / OMEGA_MIN_K17 - 0.5) / 3))
            prm["model"] = "cumulant"
        else:
            prm.update(k17_parameters(omega))
    if prm["model"] in ("cumulant", "bgk"):
        prm.update(w2=1.0, w3=1.0, w4=1.0, w5=1.0, A=0.0, B=0.0)
    if prm["model"] == "bgk":
        prm["w2"] = omega
    return prm


def collide(f, force, prm):
    """Colisao in-place.  Devolve (rho, ux, uy, uz) pre-colisao.

    Passos [G15, Sec. 4 e Apendice J; G17, Apendice]:
      1. rho, u = (sum f c + F/2)/rho
      2. f -> momentos centrais k (normalizados por rho)
      3. cumulantes  c_2 = k_2, c_3 = k_3,
         c_4: c211 = k211 - (k200 k011 + 2 k110 k101), c220 = k220 - (k200 k020 + 2 k110^2)
         c_5, c_6: ver back_transform (relaxam para 0)
      4. relaxacao (ver abaixo), 5. cumulantes -> momentos centrais -> f
    """
    w1 = prm["w1"]
    rho, ux, uy, uz = macroscopic(f, force)
    _fwd(f, uz, 2)
    _fwd(f, uy, 1)
    _fwd(f, ux, 0)
    f /= rho
    k = f

    if prm["model"] == "bgk":
        # BGK = relaxacao de todos os momentos centrais com a mesma taxa
        k1x, k1y, k1z = -k[1, 0, 0], -k[0, 1, 0], -k[0, 0, 1]
        k += w1 * (KEQ.reshape((3, 3, 3) + (1,) * rho.ndim) - k)
        k[1, 0, 0], k[0, 1, 0], k[0, 0, 1] = k1x, k1y, k1z
    else:
        _cumulant_relax(k, ux, uy, uz, prm)

    f *= rho
    _bwd(f, ux, 0)
    _bwd(f, uy, 1)
    _bwd(f, uz, 2)
    return rho, ux, uy, uz


def _cumulant_relax(k, ux, uy, uz, prm):
    w1, w2, w3, w4, w5 = prm["w1"], prm["w2"], prm["w3"], prm["w4"], prm["w5"]
    lim = prm["limiter"]

    k200, k020, k002 = k[2, 0, 0], k[0, 2, 0], k[0, 0, 2]
    k110, k101, k011 = k[1, 1, 0], k[1, 0, 1], k[0, 1, 1]

    # --- gradientes de velocidade a partir dos cumulantes de 2a ordem [G15, Eqs. 58-60]
    T = k200 + k020 + k002
    dxux = -0.5 * w1 * (2.0 * k200 - k020 - k002) - 0.5 * w2 * (T - 1.0)
    dyuy = dxux + 1.5 * w1 * (k200 - k020)
    dzuz = dxux + 1.5 * w1 * (k200 - k002)
    dxy = -3.0 * w1 * k110
    dxz = -3.0 * w1 * k101
    dyz = -3.0 * w1 * k011

    # --- 2a ordem: partes desviadoras (w1) e traco (w2) + correcao galileana [G15, Eqs. 61-63]
    ux2, uy2, uz2 = ux * ux, uy * uy, uz * uz
    mxxMyy = (1.0 - w1) * (k200 - k020) - 3.0 * (1.0 - 0.5 * w1) * (ux2 * dxux - uy2 * dyuy)
    mxxMzz = (1.0 - w1) * (k200 - k002) - 3.0 * (1.0 - 0.5 * w1) * (ux2 * dxux - uz2 * dzuz)
    T = T + w2 * (1.0 - T) - 3.0 * (1.0 - 0.5 * w2) * (ux2 * dxux + uy2 * dyuy + uz2 * dzuz)
    p200 = (T + mxxMyy + mxxMzz) / 3.0
    p020 = (T - 2.0 * mxxMyy + mxxMzz) / 3.0
    p002 = (T + mxxMyy - 2.0 * mxxMzz) / 3.0
    p110 = (1.0 - w1) * k110
    p101 = (1.0 - w1) * k101
    p011 = (1.0 - w1) * k011

    # --- 3a ordem: combinacoes simetricas (w3), anti-simetricas (w4) e c111 (w5) [G17]
    def relax(x, w):
        if lim is not None and w != 1.0:
            ax = np.abs(x)
            w = w + (1.0 - w) * ax / (lim + ax)  # limitador opcional [G17, Eq. 20]
        return (1.0 - w) * x

    def pair(a, b):
        P = relax(a + b, w3)
        M = relax(a - b, w4)
        return 0.5 * (P + M), 0.5 * (P - M)

    p210, p012 = pair(k[2, 1, 0], k[0, 1, 2])
    p201, p021 = pair(k[2, 0, 1], k[0, 2, 1])
    p120, p102 = pair(k[1, 2, 0], k[1, 0, 2])
    p111 = relax(k[1, 1, 1], w5)

    # --- 4a ordem: equilibrio parametrizado [G17] (taxa = 1); 0 no cumulante classico
    if prm["model"] == "k17":
        fac = -(1.0 / w1 - 0.5)
        A, B = prm["A"], prm["B"]
        c022 = fac * (2.0 / 3.0) * A * (dyuy + dzuz)
        c202 = fac * (2.0 / 3.0) * A * (dxux + dzuz)
        c220 = fac * (2.0 / 3.0) * A * (dxux + dyuy)
        c112 = fac * (1.0 / 3.0) * B * dxy
        c121 = fac * (1.0 / 3.0) * B * dxz
        c211 = fac * (1.0 / 3.0) * B * dyz
    else:
        c022 = c202 = c220 = c112 = c121 = c211 = 0.0

    # --- cumulantes pos-colisao -> momentos centrais (5a e 6a ordem relaxam para 0)
    p211 = c211 + p200 * p011 + 2.0 * p110 * p101
    p121 = c121 + p020 * p101 + 2.0 * p110 * p011
    p112 = c112 + p002 * p110 + 2.0 * p101 * p011
    p220 = c220 + p200 * p020 + 2.0 * p110 * p110
    p202 = c202 + p200 * p002 + 2.0 * p101 * p101
    p022 = c022 + p020 * p002 + 2.0 * p011 * p011

    p221 = (p200 * p021 + p020 * p201 + 4.0 * p110 * p111
            + 2.0 * p101 * p120 + 2.0 * p011 * p210)
    p212 = (p200 * p012 + p002 * p210 + 4.0 * p101 * p111
            + 2.0 * p110 * p102 + 2.0 * p011 * p201)
    p122 = (p002 * p120 + p020 * p102 + 4.0 * p011 * p111
            + 2.0 * p101 * p021 + 2.0 * p110 * p012)

    S4 = (p200 * p022 + p020 * p202 + p002 * p220
          + 4.0 * (p110 * p112 + p101 * p121 + p011 * p211))
    S33 = 2.0 * (p210 * p012 + p201 * p021 + p120 * p102) + 4.0 * p111 * p111
    S222 = (p200 * p020 * p002
            + 2.0 * (p200 * p011 * p011 + p020 * p101 * p101 + p002 * p110 * p110)
            + 8.0 * p110 * p101 * p011)
    p222 = S4 + S33 - 2.0 * S222

    # forca de corpo: momento central de 1a ordem troca de sinal (-F/2rho -> +F/2rho) [G15]
    k[1, 0, 0] = -k[1, 0, 0]
    k[0, 1, 0] = -k[0, 1, 0]
    k[0, 0, 1] = -k[0, 0, 1]
    k[0, 0, 0] = 1.0
    k[2, 0, 0], k[0, 2, 0], k[0, 0, 2] = p200, p020, p002
    k[1, 1, 0], k[1, 0, 1], k[0, 1, 1] = p110, p101, p011
    k[2, 1, 0], k[0, 1, 2], k[2, 0, 1], k[0, 2, 1] = p210, p012, p201, p021
    k[1, 2, 0], k[1, 0, 2], k[1, 1, 1] = p120, p102, p111
    k[2, 2, 0], k[2, 0, 2], k[0, 2, 2] = p220, p202, p022
    k[2, 1, 1], k[1, 2, 1], k[1, 1, 2] = p211, p121, p112
    k[2, 2, 1], k[2, 1, 2], k[1, 2, 2] = p221, p212, p122
    k[2, 2, 2] = p222


def stream(src, dst):
    for a in range(3):
        for b in range(3):
            for c in range(3):
                dst[a, b, c] = np.roll(src[a, b, c], (a - 1, b - 1, c - 1), axis=(0, 1, 2))


# ---------------------------------------------------------------------------
# Condicoes de contorno de parede (lista de "links" fluido -> solido)
#   f_dbar(xf, t+1) = A1 f*_d(xf) + A2 f*_d(xf - c_d) + A3 f*_dbar(xf) + mw
#   q = 1/2 (meio caminho) : A1 = 1                                  [L94]
#   q < 1/2                : A1 = 2q, A2 = 1 - 2q                    [B01]
#   q >= 1/2               : A1 = 1/(2q), A3 = (2q - 1)/(2q)         [B01]
#   parede movel           : mw = -6 w_d rho_w (c_d . u_w)           [L94]
# ---------------------------------------------------------------------------
class Links:
    def __init__(self):
        self.n = np.zeros(0, np.int64)
        self.d = np.zeros(0, np.int64)
        self.n2 = np.zeros(0, np.int64)
        self.A1 = np.zeros(0)
        self.A2 = np.zeros(0)
        self.A3 = np.zeros(0)
        self.mw = np.zeros(0)
        self.tag = np.zeros(0, np.int8)

    def apply(self, fpost, fnew):
        if self.n.size == 0:
            return
        P = fpost.reshape(27, -1)
        N = fnew.reshape(27, -1)
        db = OPP[self.d]
        N[db, self.n] = (self.A1 * P[self.d, self.n] + self.A2 * P[self.d, self.n2]
                         + self.A3 * P[db, self.n] + self.mw)

    def force(self, fpost, fnew, tag):
        """Forca no corpo (troca de quantidade de movimento) [M02]."""
        sel = self.tag == tag
        P = fpost.reshape(27, -1)
        N = fnew.reshape(27, -1)
        d, n = self.d[sel], self.n[sel]
        g = P[d, n] + N[OPP[d], n]
        return (CV[d] * g[:, None]).sum(axis=0)


def build_links(solid, wall_velocity=None, curved=None):
    """solid: mascara booleana.  wall_velocity(i,j,k,xm,ym,zm) -> (ux,uy,uz) no
    ponto medio do link.  curved: dict(mask=..., q=func(xf,yf,zf,cx,cy,cz), tag=1)."""
    L = Links()
    shape = solid.shape
    fluid = ~solid
    idx = np.arange(solid.size).reshape(shape)
    lists = {k: [] for k in ("n", "d", "n2", "A1", "A2", "A3", "mw", "tag")}
    for d in range(27):
        if d == 13:
            continue
        c = CV[d]
        nb_solid = np.roll(solid, tuple(-c), axis=(0, 1, 2))  # solid[x + c]
        m = fluid & nb_solid
        if not m.any():
            continue
        ii, jj, kk = np.nonzero(m)
        si, sj, sk = (ii + c[0]) % shape[0], (jj + c[1]) % shape[1], (kk + c[2]) % shape[2]
        bi, bj, bk = (ii - c[0]) % shape[0], (jj - c[1]) % shape[1], (kk - c[2]) % shape[2]
        n = idx[ii, jj, kk]
        n2 = idx[bi, bj, bk]
        back_fluid = fluid[bi, bj, bk]
        A1 = np.ones(n.size)
        A2 = np.zeros(n.size)
        A3 = np.zeros(n.size)
        tag = np.zeros(n.size, np.int8)
        if curved is not None:
            cm = curved["mask"][si, sj, sk]
            if cm.any():
                q = curved["q"](ii[cm] - 0.5, jj[cm] - 0.5, kk[cm] - 0.5, *c)
                a1 = np.where(q < 0.5, 2 * q, 1.0 / (2 * q))
                a2 = np.where(q < 0.5, 1 - 2 * q, 0.0)
                a3 = np.where(q < 0.5, 0.0, (2 * q - 1) / (2 * q))
                bad = (q < 0.5) & ~back_fluid[cm]  # vizinho atras nao e fluido
                a1[bad], a2[bad] = 1.0, 0.0
                A1[cm], A2[cm], A3[cm] = a1, a2, a3
                tag[cm] = curved.get("tag", 1)
        mw = np.zeros(n.size)
        if wall_velocity is not None:
            uw = wall_velocity(si, sj, sk, ii - 0.5 + 0.5 * c[0], jj - 0.5 + 0.5 * c[1],
                               kk - 0.5 + 0.5 * c[2])
            mw = -6.0 * W.ravel()[d] * (c[0] * uw[0] + c[1] * uw[1] + c[2] * uw[2])
            mw = np.broadcast_to(mw, n.shape).copy()
        for key, val in (("n", n), ("d", np.full(n.size, d)), ("n2", n2), ("A1", A1),
                         ("A2", A2), ("A3", A3), ("mw", mw), ("tag", tag)):
            lists[key].append(val)
    for key, val in lists.items():
        if val:
            setattr(L, key, np.concatenate(val))
    return L


# ---------------------------------------------------------------------------
# Casos
# ---------------------------------------------------------------------------
class Case:
    """Classe base.  Coordenadas de rede: no (i,j,k) fica em (i-0.5, j-0.5, k-0.5)
    quando ha camada fantasma (paredes no meio do caminho entre nos)."""
    name = "base"
    steady = True
    slice_axis = 2

    def __init__(self, args):
        self.args = args
        self.force = (0.0, 0.0, 0.0)

    def post_stream(self, sim, fpost, fnew):
        pass

    def report(self, sim):
        return {}

    def ref_ranges(self):
        """Faixas de referencia {grandeza: (min, max)} mostradas no painel HTML."""
        return {}

    def finish(self, sim, outdir):
        return ""


class ShearWave(Case):
    """Onda de cisalhamento: u_x = U0 sin(2 pi y / N), decai com exp(-nu k^2 t)."""
    name = "shearwave"
    steady = False
    slice_axis = 0

    def __init__(self, args):
        super().__init__(args)
        N = args.res or 32
        self.N = N
        self.nu = args.nu if args.nu else 0.01
        self.omega = 1.0 / (3.0 * self.nu + 0.5)
        self.shape = (4, N, 4)
        self.solid = np.zeros(self.shape, bool)
        self.links = Links()
        self.kw = 2 * np.pi / N
        self.U0 = 0.01
        self.ushift = args.u_adv  # velocidade de adveccao (teste de invariancia galileana)
        y = np.arange(N)
        self.sy = np.sin(self.kw * y)
        self.decay_time = 1.0 / (self.nu * self.kw ** 2)
        if not args.max_iter:
            args.max_iter = int(round(self.decay_time / 100.0)) * 100

    def initial(self):
        rho = np.ones(self.shape)
        ux = self.U0 * np.broadcast_to(self.sy[None, :, None], self.shape).copy()
        uz = np.full(self.shape, self.ushift)
        return rho, ux, np.zeros(self.shape), uz

    def report(self, sim):
        ux = sim.u[0].mean(axis=(0, 2))
        return {"amplitude": 2.0 * np.mean(ux * self.sy)}

    def finish(self, sim, outdir):
        t = np.array(sim.hist["it"], float)
        A = np.array(sim.hist["amplitude"])
        sel = t >= 0.1 * self.decay_time
        if sel.sum() < 2:
            return "Poucos pontos para medir a viscosidade."
        slope = np.polyfit(t[sel], np.log(A[sel]), 1)[0]
        nu_num = -slope / self.kw ** 2
        err = (nu_num - self.nu) / self.nu
        return ("Onda de cisalhamento N=%d, modelo=%s: nu_teorico=%.6e  nu_medido=%.6e  "
                "erro relativo=%.3e" % (self.N, sim.prm["model"], self.nu, nu_num, err))


class Duct(Case):
    """Duto retangular periodico em x, forcado por forca de corpo (gradiente de pressao)."""
    name = "duct"
    slice_axis = 0

    def __init__(self, args):
        super().__init__(args)
        Ny = args.res or 24
        Nz = int(round(Ny * 2 / 3)) if not args.res_z else args.res_z
        self.Ny, self.Nz = Ny, Nz
        self.nu = args.nu if args.nu else 0.02
        self.omega = 1.0 / (3.0 * self.nu + 0.5)
        self.shape = (4, Ny + 2, Nz + 2)
        s = np.zeros(self.shape, bool)
        s[:, 0, :] = s[:, -1, :] = s[:, :, 0] = s[:, :, -1] = True
        self.solid = s
        self.links = build_links(s)
        # solucao analitica com G = 1, depois escalada para u_max = 0.05
        self.a, self.b = Ny / 2.0, Nz / 2.0
        y = np.arange(Ny + 2) - 0.5 - self.a
        z = np.arange(Nz + 2) - 0.5 - self.b
        self.Y, self.Z = np.meshgrid(y, z, indexing="ij")
        u1 = self.analytic(1.0)
        G = 0.05 / u1.max()
        self.G = G
        self.force = (G, 0.0, 0.0)
        self.uex = self.analytic(G)

    def analytic(self, G, nterms=101):
        """F.M. White, Viscous Fluid Flow (2006), Eq. (3-48), -a<y<a, -b<z<b."""
        a, b, mu = self.a, self.b, self.nu
        u = np.zeros_like(self.Y)
        for i in range(1, 2 * nterms, 2):
            u += ((-1) ** ((i - 1) // 2) * (1.0 - np.cosh(i * np.pi * self.Z / (2 * a))
                                             / np.cosh(i * np.pi * b / (2 * a)))
                  * np.cos(i * np.pi * self.Y / (2 * a)) / i ** 3)
        u *= 16.0 * a * a * G / (mu * np.pi ** 3)
        u[0, :] = u[-1, :] = u[:, 0] = u[:, -1] = 0.0
        return u

    def initial(self):
        z = np.zeros(self.shape)
        return np.ones(self.shape), z, z.copy(), z.copy()

    def error(self, sim):
        ux = sim.u[0].mean(axis=0)[1:-1, 1:-1]
        ue = self.uex[1:-1, 1:-1]
        return np.sqrt(np.sum((ux - ue) ** 2) / np.sum(ue ** 2)), ux, ue

    def report(self, sim):
        return {"L2_error": self.error(sim)[0]}

    def finish(self, sim, outdir):
        err, ux, ue = self.error(sim)
        try:
            import matplotlib.pyplot as plt
            fig, ax = plt.subplots(1, 2, figsize=(10, 4))
            j = self.Ny // 2
            yy = np.arange(1, self.Ny + 1) - 0.5 - self.a
            zz = np.arange(1, self.Nz + 1) - 0.5 - self.b
            ax[0].plot(yy, ue[:, self.Nz // 2], "k-", label="Analitico (White 2006)")
            ax[0].plot(yy, ux[:, self.Nz // 2], "ro", mfc="none", label="LBM")
            ax[0].set_xlabel("y"); ax[0].set_ylabel("u_x"); ax[0].legend()
            ax[1].plot(zz, ue[j, :], "k-", label="Analitico")
            ax[1].plot(zz, ux[j, :], "bs", mfc="none", label="LBM")
            ax[1].set_xlabel("z"); ax[1].legend()
            fig.suptitle("Duto retangular %dx%d - erro L2 = %.2e" % (self.Ny, self.Nz, err))
            fig.tight_layout()
            fig.savefig(os.path.join(outdir, "duct_validacao.png"), dpi=120)
            plt.close(fig)
        except Exception as e:  # pragma: no cover
            print("grafico nao gerado:", e)
        return ("Duto %dx%d, nu=%.3f, modelo=%s: erro L2 relativo = %.3e, "
                "erro max = %.3e" % (self.Ny, self.Nz, self.nu, sim.prm["model"], err,
                                     np.abs(ux - ue).max() / ue.max()))


class Cylinder(Case):
    """Schaefer & Turek (1996), caso 3D-1Z: cilindro quadrado-secao circular num canal
    de secao quadrada H=0.41 m, D=0.1 m, Re = U_media D / nu = 20 (U_media = 4 Um / 9)."""
    name = "cylinder"
    slice_axis = 2
    measure_force = True
    REF = {"C_D": (6.05, 6.25), "C_L": (0.008, 0.010), "dp": (0.165, 0.175)}

    def __init__(self, args):
        super().__init__(args)
        D = args.res or 10
        H = int(round(4.1 * D))
        self.D, self.H = D, H
        Lx = int(round(25 * D))
        self.Re = args.Re if args.Re else 20.0
        self.Um = args.u_lattice if args.u_lattice else 0.1  # velocidade maxima na entrada
        self.Ubar = 4.0 * self.Um / 9.0
        self.nu = self.Ubar * D / self.Re
        self.omega = 1.0 / (3.0 * self.nu + 0.5)
        self.shape = (Lx + 1, H + 2, H + 2)
        self.xc, self.yc, self.R = 5.0 * D, 2.0 * D, 0.5 * D
        nx, ny, nz = self.shape
        X, Y = np.meshgrid(np.arange(nx) - 0.5, np.arange(ny) - 0.5, indexing="ij")
        cyl2d = (X - self.xc) ** 2 + (Y - self.yc) ** 2 < self.R ** 2
        cyl = np.broadcast_to(cyl2d[:, :, None], self.shape).copy()
        cyl[:, :, 0] = cyl[:, :, -1] = False
        s = np.zeros(self.shape, bool)
        s[0] = True                                   # entrada (parede movel)
        s[:, 0, :] = s[:, -1, :] = s[:, :, 0] = s[:, :, -1] = True
        self.cyl = cyl
        self.solid = s | cyl

        def wall_velocity(si, sj, sk, xm, ym, zm):
            u = self.inflow(ym, zm) * ((si == 0) & ~self.cyl[si, sj, sk])
            return (u, 0.0, 0.0)

        def qfun(x, y, z, cx, cy, cz):
            px, py = x - self.xc, y - self.yc
            a = cx * cx + cy * cy
            b = px * cx + py * cy
            c = px * px + py * py - self.R ** 2
            return (-b - np.sqrt(np.maximum(b * b - a * c, 0.0))) / a

        self.links = build_links(self.solid, wall_velocity,
                                 dict(mask=cyl, q=qfun, tag=1))

    def inflow(self, y, z):
        H = self.H
        y = np.clip(y, 0, H)
        z = np.clip(z, 0, H)
        return 16.0 * self.Um * y * z * (H - y) * (H - z) / H ** 4

    def initial(self):
        nx, ny, nz = self.shape
        y = np.arange(ny) - 0.5
        z = np.arange(nz) - 0.5
        u = np.broadcast_to(self.inflow(y[:, None], z[None, :])[None], self.shape).copy()
        u[self.solid] = 0.0
        zr = np.zeros(self.shape)
        return np.ones(self.shape), u, zr, zr.copy()

    def post_stream(self, sim, fpost, fnew):
        # saida: extrapolacao de nao-equilibrio com rho = 1 [G02]
        fb = fnew[:, :, :, -2]
        rho, ux, uy, uz = macroscopic(fb)
        fnew[:, :, :, -1] = equilibrium(np.ones_like(rho), ux, uy, uz) + fb - equilibrium(rho, ux, uy, uz)

    def _sample_line(self, rho, i, y, z):
        j0, k0 = int(np.floor(y + 0.5)), int(np.floor(z + 0.5))
        wy, wz = y + 0.5 - j0, z + 0.5 - k0
        return ((1 - wy) * (1 - wz) * rho[i, j0, k0] + wy * (1 - wz) * rho[i, j0 + 1, k0]
                + (1 - wy) * wz * rho[i, j0, k0 + 1] + wy * wz * rho[i, j0 + 1, k0 + 1])

    def coefficients(self, sim):
        F = sim.body_force
        Fd, Fl = F[0], F[1]
        norm = 0.5 * self.Ubar ** 2 * self.D * self.H
        rho = sim.rho
        zc = self.H / 2.0
        xf, xb = self.xc - self.R, self.xc + self.R
        i1 = int(np.ceil(xf + 0.5)) - 1
        r1, r0 = self._sample_line(rho, i1, self.yc, zc), self._sample_line(rho, i1 - 1, self.yc, zc)
        rf = r1 + (xf - (i1 - 0.5)) * (r1 - r0)
        i1 = int(np.floor(xb + 0.5)) + 1
        r1, r2 = self._sample_line(rho, i1, self.yc, zc), self._sample_line(rho, i1 + 1, self.yc, zc)
        rb = r1 - ((i1 - 0.5) - xb) * (r2 - r1)
        dp_phys = (rf - rb) / 3.0 * (0.2 / self.Ubar) ** 2  # U_media fisica = 0.2 m/s, rho = 1
        return {"C_D": Fd / norm, "C_L": Fl / norm, "dp": dp_phys}

    def report(self, sim):
        return self.coefficients(sim)

    def ref_ranges(self):
        return self.REF

    def finish(self, sim, outdir):
        c = self.coefficients(sim)
        lines = ["Schaefer & Turek 3D-1Z (Re=%.0f), D=%d celulas, modelo=%s:"
                 % (self.Re, self.D, sim.prm["model"])]
        for key, label in (("C_D", "C_D"), ("C_L", "C_L"), ("dp", "Delta p [Pa]")):
            lo, hi = self.REF[key]
            ok = "OK" if lo <= c[key] <= hi else "fora da faixa"
            lines.append("  %-13s = %.5f   faixa de referencia [%.3f, %.3f]  -> %s"
                         % (label, c[key], lo, hi, ok))
        return "\n".join(lines)


class Cavity(Case):
    """Cavidade cubica com tampa (z = N) movendo-se em +x."""
    name = "cavity"
    slice_axis = 1

    def __init__(self, args):
        super().__init__(args)
        N = args.res or 32
        self.N = N
        self.U = args.u_lattice if args.u_lattice else 0.1
        self.Re = args.Re if args.Re else 100.0
        self.nu = self.U * N / self.Re
        self.omega = 1.0 / (3.0 * self.nu + 0.5)
        self.shape = (N + 2, N + 2, N + 2)
        s = np.ones(self.shape, bool)
        s[1:-1, 1:-1, 1:-1] = False
        self.solid = s

        def wall_velocity(si, sj, sk, xm, ym, zm):
            lid = (sk == N + 1) & (si >= 1) & (si <= N) & (sj >= 1) & (sj <= N)
            return (self.U * lid, 0.0, 0.0)

        self.links = build_links(s, wall_velocity)

    def initial(self):
        z = np.zeros(self.shape)
        return np.ones(self.shape), z, z.copy(), z.copy()

    def finish(self, sim, outdir):
        N, c = self.N, self.N // 2 + 1
        ux = sim.u[0][c, c, 1:-1] / self.U
        uz = sim.u[2][1:-1, c, c] / self.U
        s = (np.arange(1, N + 1) - 0.5) / N
        np.savetxt(os.path.join(outdir, "cavity_centerlines.txt"), np.c_[s, ux, uz],
                   header="s  u_x(x=0.5,y=0.5,z=s)/U  u_z(x=s,y=0.5,z=0.5)/U")
        return ("Cavidade Re=%.0f: perfis centrais salvos em cavity_centerlines.txt "
                "(compare com Albensoeder & Kuhlmann 2005 / Ku et al. 1987)" % self.Re)


class SphereArray(Case):
    """Escoamento de Stokes atraves de um arranjo cubico simples de esferas (celula
    periodica L^3, uma esfera de diametro d = chi L), forcado por aceleracao g em x.

    Referencia analitica: A.S. Sangani, A. Acrivos, "Slow flow through a periodic array
    of spheres", Int. J. Multiphase Flow 8 (1982) 343-360 (estende Hasimoto, J. Fluid
    Mech. 5 (1959) 317).  Arraste adimensional, com u_s = velocidade superficial
    (media sobre L^3) e f_t = forca total por esfera (arraste + empuxo V_p rho g):
        C = f_t / (3 pi mu d u_s)
    Valores de C(chi) tabelados em Bogner, Mohanty, Ruede, Int. J. Multiphase Flow 68
    (2015) 71-79, Tab. 1 (arXiv:1401.2025), recuperados de C*/(1 + erro) nas duas
    resolucoes.  Serie de Sangani-Acrivos (phi = pi chi^3 / 6, x = phi^(1/3)):
        1/C = 1 - 1.7601 x + x^3 - 1.5593 x^6 + 3.9799 x^8 - 3.0734 x^10 + ...
    """
    name = "spheres"
    slice_axis = 1
    measure_force = True
    REF = {0.2: 1.385, 0.3: 1.70, 0.4: 2.150, 0.5: 2.845, 0.6: 3.974, 0.7: 6.00,
           0.8: 10.05, 0.9: 19.16}

    def __init__(self, args):
        super().__init__(args)
        L = args.res or 32
        self.L = L
        self.chi = args.chi
        self.R = 0.5 * self.chi * L
        self.d = 2.0 * self.R
        self.nu = args.nu if args.nu else 0.02
        self.omega = 1.0 / (3.0 * self.nu + 0.5)
        self.shape = (L, L, L)
        x = np.arange(L) - 0.5
        X, Y, Z = np.meshgrid(x, x, x, indexing="ij")
        self.c0 = 0.5 * L - 0.5  # centro entre nos (geometria simetrica)
        sph = (X - self.c0) ** 2 + (Y - self.c0) ** 2 + (Z - self.c0) ** 2 < self.R ** 2
        self.solid = sph
        self.phi = np.pi * self.chi ** 3 / 6.0
        self.Vp = np.pi * self.d ** 3 / 6.0

        def qfun(x, y, z, cx, cy, cz):
            px, py, pz = x - self.c0, y - self.c0, z - self.c0
            a = cx * cx + cy * cy + cz * cz
            b = px * cx + py * cy + pz * cz
            c = px * px + py * py + pz * pz - self.R ** 2
            return (-b - np.sqrt(np.maximum(b * b - a * c, 0.0))) / a

        self.links = build_links(sph, curved=dict(mask=sph, q=qfun, tag=1))
        # aceleracao escolhida para u_s ~ 1e-4 (Re_p = u_s d / nu << 1, regime de Stokes)
        Cref = self.reference()
        self.g = 1e-4 * 3.0 * np.pi * self.nu * self.d * Cref / L ** 3
        self.force = (self.g, 0.0, 0.0)

    def reference(self):
        if round(self.chi, 3) in self.REF:
            return self.REF[round(self.chi, 3)]
        return self.series()

    def ref_ranges(self):
        r = self.reference()
        return {"C": (r, r), "C_balanco": (r, r)}

    def series(self):
        x = self.phi ** (1.0 / 3.0)
        return 1.0 / (1 - 1.7601 * x + x ** 3 - 1.5593 * x ** 6 + 3.9799 * x ** 8
                      - 3.0734 * x ** 10)

    def initial(self):
        z = np.zeros(self.shape)
        return np.ones(self.shape), z, z.copy(), z.copy()

    def drag(self, sim):
        us = sim.u[0][sim.fluid].sum() / self.L ** 3
        ft_me = sim.body_force[0] + self.g * self.Vp      # troca de q. de movimento + empuxo
        ft_bal = self.g * self.L ** 3                      # balanco global (estado estacionario)
        norm = 3.0 * np.pi * self.nu * self.d * us
        return ft_me / norm, ft_bal / norm, us

    def report(self, sim):
        C, Cbal, us = self.drag(sim)
        return {"C": C, "C_balanco": Cbal, "Re_p": us * self.d / self.nu}

    def finish(self, sim, outdir):
        C, Cbal, us = self.drag(sim)
        ref = self.reference()
        return ("Arranjo cubico de esferas (Sangani & Acrivos 1982), L=%d, chi=%.2f "
                "(d=%.1f, phi=%.4f), modelo=%s:\n"
                "  C (troca de q. mov.) = %.4f   C (balanco global) = %.4f\n"
                "  referencia C = %.4f  -> erro = %+.2f%%   (serie S&A: %.4f)\n"
                "  Re_p = %.2e,  fracao solida discreta = %.4f"
                % (self.L, self.chi, self.d, self.phi, sim.prm["model"], C, Cbal, ref,
                   100 * (C - ref) / ref, self.series(), us * self.d / self.nu,
                   self.solid.mean()))


CASES = {"shearwave": ShearWave, "duct": Duct, "cylinder": Cylinder, "cavity": Cavity,
         "spheres": SphereArray}


# ---------------------------------------------------------------------------
# Simulacao
# ---------------------------------------------------------------------------
class Simulation:
    def __init__(self, case, args):
        self.case = case
        self.args = args
        self.prm = collision_parameters(case.omega, args.model, args.limiter)
        self.solid_idx = np.flatnonzero(case.solid)
        self.fluid = ~case.solid
        rho, ux, uy, uz = case.initial()
        self.f = equilibrium(rho, ux, uy, uz)
        self.ftmp = np.empty_like(self.f)
        self.it = 0
        self.hist = {"it": [], "res": []}
        self.rho, self.u = rho, np.stack([ux, uy, uz])
        self.u_prev = self.u.copy()
        self.body_force = np.zeros(3)

    def step(self):
        case = self.case
        # colisao em fatias de planos x: os temporarios do numpy cabem no cache (~3x mais rapido)
        nx, ny, nz = case.shape
        ch = max(1, 4096 // (ny * nz))
        for x0 in range(0, nx, ch):
            r, ux, uy, uz = collide(self.f[..., x0:x0 + ch, :, :], case.force, self.prm)
            self.rho[x0:x0 + ch] = r
            self.u[0, x0:x0 + ch], self.u[1, x0:x0 + ch], self.u[2, x0:x0 + ch] = ux, uy, uz
        self.f.reshape(27, -1)[:, self.solid_idx] = W.reshape(27, 1)
        stream(self.f, self.ftmp)
        case.links.apply(self.f, self.ftmp)
        case.post_stream(self, self.f, self.ftmp)
        if getattr(case, "measure_force", False):
            self.body_force = case.links.force(self.f, self.ftmp, 1)
        self.f, self.ftmp = self.ftmp, self.f
        self.it += 1

    def update_macros(self):
        rho, ux, uy, uz = macroscopic(self.f, self.case.force)
        self.rho, self.u = rho, np.stack([ux, uy, uz])

    def record(self):
        """Residuo de convergencia ||u^n - u^(n-100)|| / ||u^n|| nos nos de fluido."""
        self.update_macros()
        du = (self.u - self.u_prev)[:, self.fluid]
        un = self.u[:, self.fluid]
        res = np.sqrt(np.sum(du * du) / max(np.sum(un * un), 1e-300))
        self.u_prev = self.u.copy()
        self.hist["it"].append(self.it)
        self.hist["res"].append(res)
        for key, val in self.case.report(self).items():
            self.hist.setdefault(key, []).append(float(val))
        return res

    # --- E/S ---------------------------------------------------------------
    def save_checkpoint(self, outdir):
        path = os.path.join(outdir, "checkpoint.npz")
        tmp = os.path.join(outdir, "checkpoint_tmp.npz")
        np.savez(tmp, f=self.f, it=self.it, u_prev=self.u_prev,
                 args=json.dumps(vars(self.args)),
                 hist=json.dumps({k: list(map(float, v)) for k, v in self.hist.items()}))
        os.replace(tmp, path)
        return path

    def load_checkpoint(self, path):
        data = np.load(path)
        self.f = data["f"].copy()
        self.ftmp = np.empty_like(self.f)
        self.it = int(data["it"])
        self.u_prev = data["u_prev"].copy()
        self.hist = {k: list(v) for k, v in json.loads(str(data["hist"])).items()}
        self.hist["it"] = [int(i) for i in self.hist["it"]]
        self.update_macros()

    def save_vtk(self, outdir):
        path = os.path.join(outdir, "%s_it%07d.vtk" % (self.case.name, self.it))
        self.update_macros()
        nx, ny, nz = self.rho.shape
        rho = np.where(self.fluid, self.rho, 1.0)
        u = np.where(self.fluid[None], self.u, 0.0)
        with open(path, "wb") as fh:
            fh.write(("# vtk DataFile Version 3.0\nLBM %s it=%d\nBINARY\n"
                      "DATASET STRUCTURED_POINTS\nDIMENSIONS %d %d %d\n"
                      "ORIGIN -0.5 -0.5 -0.5\nSPACING 1 1 1\nPOINT_DATA %d\n"
                      % (self.case.name, self.it, nx, ny, nz, rho.size)).encode())
            fh.write(b"SCALARS density float 1\nLOOKUP_TABLE default\n")
            fh.write(rho.transpose(2, 1, 0).astype(">f4").tobytes())
            fh.write(b"\nSCALARS solid float 1\nLOOKUP_TABLE default\n")
            fh.write((~self.fluid).transpose(2, 1, 0).astype(">f4").tobytes())
            fh.write(b"\nVECTORS velocity float\n")
            fh.write(u.transpose(3, 2, 1, 0).astype(">f4").tobytes())
            fh.write(b"\n")
        return path


# ---------------------------------------------------------------------------
# Interface grafica: curva de convergencia + campo de velocidade + controles
# ---------------------------------------------------------------------------
class Monitor:
    def __init__(self, sim, args, outdir):
        import matplotlib
        if args.no_gui:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from matplotlib.widgets import Button, CheckButtons
        self.plt = plt
        self.sim, self.outdir, self.gui = sim, outdir, not args.no_gui
        self.paused = False
        self.stop = False
        self.save_now = False
        self.autosave = args.save_every > 0
        self.fig = plt.figure(figsize=(12, 5.5))
        self.ax1 = self.fig.add_axes([0.06, 0.2, 0.42, 0.7])
        self.ax2 = self.fig.add_axes([0.55, 0.2, 0.42, 0.7])
        self.line, = self.ax1.semilogy([], [], "b.-")
        self.ax1.set_xlabel("iteracao")
        self.ax1.set_ylabel(r"$\|u^n-u^{n-100}\|_2/\|u^n\|_2$")
        self.ax1.grid(True, which="both", alpha=0.3)
        self.img = None
        self.cbar = None
        if self.gui:
            self.bp = Button(self.fig.add_axes([0.06, 0.03, 0.14, 0.07]), "Pausar (P)")
            self.bs = Button(self.fig.add_axes([0.22, 0.03, 0.14, 0.07]), "Salvar (S)")
            self.bq = Button(self.fig.add_axes([0.38, 0.03, 0.18, 0.07]), "Salvar e sair (Q)")
            self.cb = CheckButtons(self.fig.add_axes([0.6, 0.02, 0.25, 0.09]),
                                   ["Auto-salvar a cada 100 it."], [self.autosave])
            self.bp.on_clicked(lambda e: self.toggle_pause())
            self.bs.on_clicked(lambda e: self.request_save())
            self.bq.on_clicked(lambda e: self.request_stop())
            self.cb.on_clicked(lambda e: self.toggle_autosave())
            self.fig.canvas.mpl_connect("key_press_event", self.on_key)
            self.fig.canvas.mpl_connect("close_event", lambda e: self.request_stop())
            plt.ion()
            plt.show(block=False)

    def toggle_pause(self):
        self.paused = not self.paused
        if self.gui:
            self.bp.label.set_text("Continuar (P)" if self.paused else "Pausar (P)")
        print("PAUSADO" if self.paused else "CONTINUANDO", "na iteracao", self.sim.it)

    def toggle_autosave(self):
        self.autosave = not self.autosave
        print("Auto-salvar:", self.autosave)

    def request_save(self):
        self.save_now = True

    def request_stop(self):
        self.stop = True

    def on_key(self, ev):
        key = (ev.key or "").lower()
        if key == "p":
            self.toggle_pause()
        elif key == "s":
            self.request_save()
        elif key == "q":
            self.request_stop()

    def events(self, dt=0.001):
        if self.gui:
            self.fig.canvas.flush_events()
            if dt > 0.01:
                time.sleep(dt)

    def update(self):
        s = self.sim
        self.line.set_data(s.hist["it"], s.hist["res"])
        self.ax1.relim()
        self.ax1.autoscale_view()
        extra = [k for k in s.hist if k not in ("it", "res")]
        txt = "  ".join("%s=%.5g" % (k, s.hist[k][-1]) for k in extra)
        self.ax1.set_title("Convergencia (it=%d)  %s" % (s.it, txt), fontsize=9)
        ax = s.case.slice_axis
        n = s.u.shape[1 + ax] // 2
        umag = np.sqrt((s.u ** 2).sum(axis=0))
        sl = np.take(np.where(s.fluid, umag, np.nan), n, axis=ax).T
        if self.img is None:
            self.img = self.ax2.imshow(sl, origin="lower", cmap="jet", aspect="equal")
            self.cbar = self.fig.colorbar(self.img, ax=self.ax2, shrink=0.8)
            names = "xyz".replace("xyz"[ax], "")
            self.ax2.set_xlabel(names[0]); self.ax2.set_ylabel(names[1])
        else:
            self.img.set_data(sl)
            self.img.set_clim(np.nanmin(sl), np.nanmax(sl))
        self.ax2.set_title("|u| no plano %s = %d" % ("xyz"[ax], n))
        if self.gui:
            self.fig.canvas.draw_idle()
            self.plt.pause(0.001)

    def savefig(self):
        self.fig.savefig(os.path.join(self.outdir, "convergencia.png"), dpi=110)


def log_convergence(outdir, sim, args, mlups):
    """Acrescenta uma linha a <outdir>/convergencia.csv a cada relatorio (100 it.).
    Tambem regrava <outdir>/convergencia.html (painel no navegador, recarrega sozinho);
de outro terminal:  python acompanhar.py <outdir> [--servir]"""
    path = os.path.join(outdir, "convergencia.csv")
    keys = [k for k in sim.hist if k not in ("it", "res")]
    new = not os.path.exists(path)
    with open(path, "a") as fh:
        if new:
            refs = "".join(" ref_%s=%.6g:%.6g" % (k, lo, hi)
                           for k, (lo, hi) in sim.case.ref_ranges().items())
            fh.write("# caso=%s tol=%g max_iter=%d steady=%d modelo=%s%s\n"
                     % (sim.case.name, args.tol, args.max_iter, int(sim.case.steady),
                        sim.prm["model"], refs))
            fh.write(",".join(["it", "res", "tempo_unix", "mlups"] + keys) + "\n")
        fh.write(",".join(["%d" % sim.it, "%.6e" % sim.hist["res"][-1], "%.3f" % time.time(),
                           "%.4f" % mlups] + ["%.8g" % sim.hist[k][-1] for k in keys]) + "\n")
    try:  # painel HTML <outdir>/convergencia.html (abra no navegador; recarrega sozinho)
        import acompanhar
        acompanhar.gerar_html(outdir)
    except Exception as e:  # o painel nunca interrompe a simulacao
        if not getattr(log_convergence, "_avisou", False):
            print("AVISO: painel HTML nao gerado (%s)" % e)
            log_convergence._avisou = True


def parse_args(argv=None):
    p = argparse.ArgumentParser(description="LBM D3Q27 cumulante (K17) para escoamentos laminares 3D")
    p.add_argument("--case", choices=sorted(CASES), default="duct")
    p.add_argument("--model", choices=["k17", "cumulant", "bgk"], default="k17",
                   help="k17 = cumulante parametrizado [G17] (padrao); cumulant = [G15]; bgk")
    p.add_argument("--res", type=int, default=0, help="resolucao (D do cilindro, N, Ny ...)")
    p.add_argument("--res-z", type=int, default=0, help="Nz do duto")
    p.add_argument("--chi", type=float, default=0.5,
                   help="spheres: razao d/L (fracao solida phi = pi chi^3 / 6)")
    p.add_argument("--Re", type=float, default=0.0)
    p.add_argument("--nu", type=float, default=0.0, help="viscosidade de rede (duct/shearwave)")
    p.add_argument("--u-lattice", type=float, default=0.0, help="velocidade de referencia de rede")
    p.add_argument("--u-adv", type=float, default=0.0, help="adveccao na onda de cisalhamento")
    p.add_argument("--limiter", type=float, default=None, help="limitador de 3a ordem [G17] (ex.: 0.01)")
    p.add_argument("--max-iter", type=int, default=0)
    p.add_argument("--tol", type=float, default=1e-7, help="criterio de parada do residuo")
    p.add_argument("--save-every", type=int, default=100,
                   help="salvar VTK + checkpoint a cada N iteracoes (0 desliga; GUI pode alternar)")
    p.add_argument("--outdir", default="")
    p.add_argument("--restart", default="", help="arquivo checkpoint.npz para reiniciar")
    p.add_argument("--no-gui", action="store_true", help="sem janela (Ctrl+C salva e sai)")
    p.add_argument("--report-every", type=int, default=100)
    p.add_argument("--max-time", type=float, default=0.0,
                   help="limite de tempo de relogio em minutos: ao atingir, salva o checkpoint "
                        "e encerra (util na nuvem; continue com --restart)")
    return p.parse_args(argv)


def main(argv=None):
    args = parse_args(argv)
    restart_file = args.restart or None
    if restart_file:
        # parametros fisicos vem do checkpoint; controles de execucao podem ser trocados
        saved = json.loads(str(np.load(restart_file)["args"]))
        given = {a[2:].split("=")[0].replace("-", "_")
                 for a in (sys.argv[1:] if argv is None else argv) if a.startswith("--")}
        runtime = {"no_gui", "max_iter", "save_every", "tol", "outdir", "report_every", "max_time"}
        for k, v in saved.items():
            if k != "restart" and not (k in runtime and k in given):
                setattr(args, k, v)
        args.restart = restart_file

    outdir = args.outdir or "saida_%s" % args.case
    os.makedirs(outdir, exist_ok=True)
    case = CASES[args.case](args)
    if not args.max_iter:
        args.max_iter = 200000
    sim = Simulation(case, args)
    if restart_file:
        sim.load_checkpoint(restart_file)
        print("Reiniciando de %s na iteracao %d" % (restart_file, sim.it))

    n_nodes = int(np.prod(case.shape))
    print("Caso: %s | malha %s (%d nos) | modelo %s | omega = %.4f | nu = %.5f"
          % (case.name, case.shape, n_nodes, sim.prm["model"], case.omega, case.nu))
    if sim.prm["model"] == "k17":
        print("  taxas K17: w3=%.4f w4=%.4f w5=%.4f A=%.4f B=%.4f"
              % (sim.prm["w3"], sim.prm["w4"], sim.prm["w5"], sim.prm["A"], sim.prm["B"]))

    mon = Monitor(sim, args, outdir)
    if not sim.hist["it"]:
        sim.u_prev = sim.u.copy()

    def sigint(signum, frame):
        print("\nInterrompido: salvando checkpoint...")
        mon.request_stop()
    signal.signal(signal.SIGINT, sigint)
    signal.signal(signal.SIGTERM, sigint)  # ex.: cancelamento na nuvem

    every = args.report_every
    t0, it0 = time.time(), sim.it
    converged = False
    while sim.it < args.max_iter and not mon.stop:
        if mon.paused:
            if mon.save_now:
                mon.save_now = False
                print("  salvo:", sim.save_vtk(outdir), sim.save_checkpoint(outdir))
            mon.events(0.1)
            continue
        sim.step()
        if sim.it % 10 == 0:
            mon.events()
        if sim.it % every == 0:
            res = sim.record()
            mlups = n_nodes * (sim.it - it0) / max(time.time() - t0, 1e-9) / 1e6
            extra = "  ".join("%s=%.6g" % (k, v[-1]) for k, v in sim.hist.items()
                              if k not in ("it", "res"))
            print("it %7d  res %.3e  %s  [%.2f MLUPS]" % (sim.it, res, extra, mlups))
            log_convergence(outdir, sim, args, mlups)
            if not np.isfinite(res):
                print("Simulacao divergiu (NaN).")
                break
            mon.update()
            if mon.autosave and args.save_every and sim.it % max(args.save_every, every) == 0:
                sim.save_vtk(outdir)
                sim.save_checkpoint(outdir)
            if case.steady and res < args.tol:
                converged = True
                print("Convergiu: residuo %.3e < tol %.1e" % (res, args.tol))
                break
            if args.max_time and time.time() - t0 > 60.0 * args.max_time:
                print("Limite de tempo (%g min) atingido: salvando checkpoint..." % args.max_time)
                break
        if mon.save_now:
            mon.save_now = False
            print("  salvo:", sim.save_vtk(outdir), sim.save_checkpoint(outdir))

    ck = sim.save_checkpoint(outdir)
    vtk = sim.save_vtk(outdir)
    if sim.it % every:
        sim.record()
    mon.update()
    mon.savefig()
    print("Checkpoint: %s   (reinicie com: python %s --restart %s)" % (ck, sys.argv[0], ck))
    print("Campo final (ParaView): %s" % vtk)
    if converged or not case.steady or sim.it >= args.max_iter:
        msg = case.finish(sim, outdir)
        if msg and case.steady and not converged:
            msg += "\n  (ATENCAO: max-iter atingido sem convergir; resultado provisorio)"
        if msg:
            print("=" * 78 + "\nVALIDACAO\n" + msg + "\n" + "=" * 78)
            with open(os.path.join(outdir, "validacao.txt"), "a") as fh:
                fh.write(msg + "\n")
    terminou = converged or not case.steady or sim.it >= args.max_iter
    res = sim.hist["res"][-1] if sim.hist["res"] else float("nan")
    with open(os.path.join(outdir, "estado.json"), "w") as fh:  # lido pelo fluxo na nuvem
        json.dump({"caso": case.name, "it": sim.it, "residuo": res, "tol": args.tol,
                   "convergiu": converged, "terminou": bool(terminou),
                   "divergiu": not np.isfinite(res)}, fh, indent=1)
    if mon.gui and not mon.stop:
        mon.plt.ioff()
        mon.plt.show()
    return sim


if __name__ == "__main__":
    main()
