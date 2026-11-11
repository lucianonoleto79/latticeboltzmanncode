# -*- coding: utf-8 -*-
"""
figuras.py - imagens de pre- e pos-processamento dos casos do lbm3d.py

Cada caso tem uma pasta (casos/<caso>/) com as subpastas:
  dominio/     geometria do dominio (vista 3D e cortes fluido/solido)
  contorno/    mapa das condicoes de contorno (cortes, vista 3D) + descricao em texto
  malha/       malha de calculo (cortes com os nos, detalhe da parede com os links) + malha.txt
  escoamento/  campos do escoamento (|u|, componentes, pressao, linhas de corrente) + VTK
  graficos/    convergencia (residuos), grandezas monitoradas, validacao, CSV e painel HTML

As figuras usam a API orientada a objetos do matplotlib (matplotlib.figure.Figure),
sem pyplot: nao interferem na janela de acompanhamento e funcionam sem tela.
"""
import os

import numpy as np

# codigos do mapa de condicoes de contorno
FLUIDO, PAREDE, MOVEL, CURVA, SAIDA, INTERNO = 0, 1, 2, 3, 4, 5
NOMES_CC = {
    FLUIDO: "fluido",
    PAREDE: "parede fixa (bounce-back, meio caminho) [L94]",
    MOVEL: "parede movel / entrada de velocidade [L94]",
    CURVA: "parede curva (bounce-back interpolado) [B01]",
    SAIDA: "saida com pressao imposta (rho = 1) [G02]",
    INTERNO: "solido (interior, sem contato com o fluido)",
}
CORES_CC = {
    FLUIDO: "#e8f1fb",
    PAREDE: "#7a7974",
    MOVEL: "#2a78d6",
    CURVA: "#d03b3b",
    SAIDA: "#1baf7a",
    INTERNO: "#3a3a38",
}
EIXOS = "xyz"


def _figura(w, h):
    from matplotlib.figure import Figure
    return Figure(figsize=(w, h))


def _salvar(fig, caminho, dpi=120):
    fig.savefig(caminho, dpi=dpi, bbox_inches="tight")
    return caminho


def _grade_cortes(case, altura_extra=0.0):
    """Figura com os tres cortes: lado a lado, ou empilhados para dominios alongados."""
    sh = case.shape
    if min(sh) == 1:           # malha 2D: um unico corte
        fig = _figura(8, 7 + altura_extra)
        return fig, [fig.add_subplot(1, 1, 1)]
    if max(sh) > 3 * min(sh):
        fig = _figura(12, 10 + altura_extra)
        return fig, [fig.add_subplot(3, 1, n + 1) for n in range(3)]
    fig = _figura(14, 4.6 + altura_extra)
    return fig, [fig.add_subplot(1, 3, n + 1) for n in range(3)]


def _cmap_cc():
    from matplotlib.colors import BoundaryNorm, ListedColormap
    cores = [CORES_CC[k] for k in sorted(CORES_CC)]
    return ListedColormap(cores), BoundaryNorm(np.arange(-0.5, len(cores) + 0.5), len(cores))


# ---------------------------------------------------------------------------
# Utilitarios de geometria
# ---------------------------------------------------------------------------
def planos(case):
    """Os tres cortes medios (eixo normal, indice).  O corte principal (slice_axis)
    vem primeiro.  Para casos com obstaculo o indice passa pelo seu centro."""
    centro = getattr(case, "centro_obstaculo", None)
    out = []
    eixos = [case.slice_axis] + [a for a in range(3) if a != case.slice_axis]
    if min(case.shape) == 1:   # malha 2D (uma camada): so o plano que contem o escoamento
        eixos = [int(np.argmin(case.shape))]
    for ax in eixos:
        n = case.shape[ax]
        idx = n // 2
        if centro is not None and centro[ax] is not None:
            idx = int(np.clip(np.floor(centro[ax] + 0.5), 0, n - 1))
        out.append((ax, idx))
    return out


def _corte(A, ax, idx):
    """Corte 2D de A (indices [i,j,k]) no plano normal a `ax`; devolve (B, eixo_h, eixo_v)
    com B[h, v]."""
    B = np.take(A, idx, axis=ax)
    h, v = [a for a in range(3) if a != ax]
    return B, h, v


def _extent(case, h, v):
    """Limites (em unidades de rede, ou fisicas se case.dx_fisico) do corte.
    Convencao do lbm3d: o no i fica em x = i - 0.5, ocupando a celula [i-1, i]."""
    dx = getattr(case, "dx_fisico", None) or 1.0
    return (-1.0 * dx, (case.shape[h] - 1) * dx, -1.0 * dx, (case.shape[v] - 1) * dx)


def _rotulo(case, a):
    un = getattr(case, "unidade", None)
    return "%s [%s]" % (EIXOS[a], un) if un else "%s [nos]" % EIXOS[a]


def mapa_contorno(case):
    """Mapa inteiro (shape da malha) com o tipo de cada no (codigos acima), obtido dos
    links fluido->solido montados em build_links (o mesmo dado que a simulacao usa)."""
    shape = case.shape
    cod = np.where(case.solid, INTERNO, FLUIDO).astype(np.int8)
    L = case.links
    if L.n.size:
        from lbm3d import CV
        ii = np.unravel_index(L.n, shape)
        c = CV[L.d]
        s = np.ravel_multi_index(tuple((ii[a] + c[:, a]) % shape[a] for a in range(3)), shape)
        tipo = np.where(L.tag == 1, CURVA, np.where(L.mw != 0, MOVEL, PAREDE))
        flat = cod.reshape(-1)
        for t in (PAREDE, MOVEL, CURVA):  # prioridade crescente
            flat[s[tipo == t]] = t
    saida = getattr(case, "mascara_saida", None)
    if saida is not None:
        cod[saida()] = SAIDA
    return cod


def q_links(case):
    """Fracao q do link ate a parede, recuperada dos coeficientes de Bouzidi [B01]."""
    L = case.links
    q = np.full(L.n.size, 0.5)
    curva = L.tag == 1
    a1 = L.A1[curva]
    q[curva] = np.where(L.A2[curva] > 0, 0.5 * a1, 0.5 / np.maximum(a1, 1e-300))
    q[curva & (L.A2 == 0) & (L.A3 == 0)] = 0.5  # recaida no meio caminho
    return q


# ---------------------------------------------------------------------------
# Dominio
# ---------------------------------------------------------------------------
def _caixa(ax3, case):
    dx = getattr(case, "dx_fisico", None) or 1.0
    lo = -1.0 * dx
    hi = [(n - 1) * dx for n in case.shape]
    for a in range(3):
        b, c = [k for k in range(3) if k != a]
        for vb in (lo, hi[b]):
            for vc in (lo, hi[c]):
                p0, p1 = [0.0] * 3, [0.0] * 3
                p0[a], p1[a] = lo, hi[a]
                p0[b] = p1[b] = vb
                p0[c] = p1[c] = vc
                ax3.plot(*zip(p0, p1), color="0.4", lw=0.8)


def _nuvem(ax3, case, cod, maxpts=9000, seed=0):
    rng = np.random.default_rng(seed)
    dx = getattr(case, "dx_fisico", None) or 1.0
    estilos = {PAREDE: 0.04, MOVEL: 0.25, SAIDA: 0.25, CURVA: 0.9}
    for t, alpha in estilos.items():
        p = np.argwhere(cod == t)
        if not len(p):
            continue
        if len(p) > maxpts:
            p = p[rng.choice(len(p), maxpts, replace=False)]
        x = (p - 0.5) * dx
        ax3.scatter(x[:, 0], x[:, 1], x[:, 2], s=2 if t != CURVA else 3, c=CORES_CC[t],
                    alpha=alpha, depthshade=False, label=NOMES_CC[t].split(" (")[0].split(" [")[0])
    ax3.set_box_aspect([max(n, 1) for n in case.shape])
    for a, f in enumerate((ax3.set_xlabel, ax3.set_ylabel, ax3.set_zlabel)):
        f(_rotulo(case, a))
    ax3.legend(loc="upper left", fontsize=7, markerscale=4)


def figuras_dominio(case, pasta):
    cod = mapa_contorno(case)
    arquivos = []
    fig = _figura(9, 6.5)
    ax3 = fig.add_subplot(projection="3d")
    _caixa(ax3, case)
    _nuvem(ax3, case, cod)
    ax3.set_title("Dominio computacional - %s  (malha %d x %d x %d)" % ((case.name,) + tuple(case.shape)))
    arquivos.append(_salvar(fig, os.path.join(pasta, "dominio_3d.png")))

    from matplotlib.colors import ListedColormap
    fig, eixos = _grade_cortes(case)
    for a, (ax, idx) in zip(eixos, planos(case)):
        B, h, v = _corte(case.solid, ax, idx)
        a.imshow(B.T.astype(float), origin="lower", extent=_extent(case, h, v),
                 cmap=ListedColormap([CORES_CC[FLUIDO], CORES_CC[INTERNO]]), vmin=0, vmax=1,
                 interpolation="nearest", aspect="equal")
        a.set_xlabel(_rotulo(case, h))
        a.set_ylabel(_rotulo(case, v))
        a.set_title("corte %s = no %d" % (EIXOS[ax], idx), fontsize=9)
    fig.suptitle("Dominio - %s: fluido (claro) e solido (escuro)" % case.name)
    arquivos.append(_salvar(fig, os.path.join(pasta, "dominio_cortes.png")))
    return arquivos


# ---------------------------------------------------------------------------
# Condicoes de contorno
# ---------------------------------------------------------------------------
def figuras_contorno(case, pasta):
    from matplotlib.patches import Patch
    cod = mapa_contorno(case)
    cmap, norm = _cmap_cc()
    presentes = sorted(set(np.unique(cod)))
    arquivos = []

    fig, eixos = _grade_cortes(case, 0.6)
    for a, (ax, idx) in zip(eixos, planos(case)):
        B, h, v = _corte(cod, ax, idx)
        a.imshow(B.T, origin="lower", extent=_extent(case, h, v), cmap=cmap, norm=norm,
                 interpolation="nearest", aspect="equal")
        a.set_xlabel(_rotulo(case, h))
        a.set_ylabel(_rotulo(case, v))
        a.set_title("corte %s = no %d" % (EIXOS[ax], idx), fontsize=9)
    per = getattr(case, "periodico", "")
    fig.legend(handles=[Patch(color=CORES_CC[t], label=NOMES_CC[t]) for t in presentes],
               loc="lower center", ncol=3, fontsize=8, frameon=False)
    fig.suptitle("Condicoes de contorno - %s%s" % (case.name, "  |  periodico em " + per if per else ""))
    vertical = max(case.shape) > 3 * min(case.shape)
    fig.subplots_adjust(bottom=0.12 if vertical else 0.22, hspace=0.45)
    arquivos.append(_salvar(fig, os.path.join(pasta, "contorno_cortes.png")))

    fig = _figura(9, 6.5)
    ax3 = fig.add_subplot(projection="3d")
    _caixa(ax3, case)
    _nuvem(ax3, case, cod)
    ax3.set_title("Condicoes de contorno - %s" % case.name)
    arquivos.append(_salvar(fig, os.path.join(pasta, "contorno_3d.png")))

    extra = getattr(case, "figura_contorno", None)
    if extra is not None:
        arquivos += extra(pasta, _figura, _salvar)

    linhas = ["Condicoes de contorno - caso %s" % case.name, "=" * 60, ""]
    linhas += case.descricao_contorno()
    linhas += ["", "Nos por tipo (mapa em contorno_cortes.png / contorno_3d.png):"]
    for t in presentes:
        linhas.append("  %-55s %10d" % (NOMES_CC[t], int((cod == t).sum())))
    L = case.links
    linhas += ["", "Links fluido->solido: %d  (paredes curvas: %d, paredes moveis: %d)"
               % (L.n.size, int((L.tag == 1).sum()), int((L.mw != 0).sum()))]
    caminho = os.path.join(pasta, "condicoes_contorno.txt")
    with open(caminho, "w", encoding="utf-8") as fh:
        fh.write("\n".join(linhas) + "\n")
    arquivos.append(caminho)
    return arquivos


# ---------------------------------------------------------------------------
# Malha
# ---------------------------------------------------------------------------
def _janela_detalhe(case, cod, ax, idx, h, v, lado_max=40):
    B, _, _ = _corte(cod, ax, idx)
    p = np.argwhere(B == CURVA)
    if len(p):
        lo, hi = p.min(axis=0) - 3, p.max(axis=0) + 4
        c = 0.5 * (lo + hi)
        meia = min(0.5 * max(hi - lo), lado_max / 2)
        lo, hi = np.floor(c - meia).astype(int), np.ceil(c + meia).astype(int)
    else:  # sem parede curva: canto do dominio
        lo = np.array([0, 0])
        hi = np.array([min(14, B.shape[0]), min(14, B.shape[1])])
    lo = np.maximum(lo, 0)
    hi = np.minimum(hi, B.shape)
    return lo, hi


def figuras_malha(case, pasta, prm=None):
    from matplotlib.patches import Patch
    cod = mapa_contorno(case)
    cmap, norm = _cmap_cc()
    arquivos = []
    ax, idx = planos(case)[0]
    B, h, v = _corte(cod, ax, idx)
    dx = getattr(case, "dx_fisico", None) or 1.0

    # 1) corte principal com as linhas da malha (cada celula = um no)
    fig = _figura(12, 5.5)
    a = fig.add_subplot(1, 1, 1)
    xs = (np.arange(B.shape[0] + 1) - 1.0) * dx
    ys = (np.arange(B.shape[1] + 1) - 1.0) * dx
    a.pcolormesh(xs, ys, B.T, cmap=cmap, norm=norm, edgecolors="0.55", linewidth=0.15)
    a.set_aspect("equal")
    a.set_xlabel(_rotulo(case, h))
    a.set_ylabel(_rotulo(case, v))
    a.set_title("Malha D3Q27 - corte %s = no %d  (%d x %d nos neste plano; malha total %d x %d x %d)"
                % ((EIXOS[ax], idx, B.shape[0], B.shape[1]) + tuple(case.shape)), fontsize=9)
    a.legend(handles=[Patch(color=CORES_CC[t], label=NOMES_CC[t]) for t in sorted(set(np.unique(B)))],
             loc="upper center", bbox_to_anchor=(0.5, -0.12), ncol=3, fontsize=7, frameon=False)
    arquivos.append(_salvar(fig, os.path.join(pasta, "malha_corte.png")))

    # 2) detalhe junto a parede: nos e links ate a parede (fracao q)
    lo, hi = _janela_detalhe(case, cod, ax, idx, h, v)
    fig = _figura(7.5, 7.5)
    a = fig.add_subplot(1, 1, 1)
    sub = B[lo[0]:hi[0], lo[1]:hi[1]]
    xs = (np.arange(lo[0], hi[0] + 1) - 1.0) * dx
    ys = (np.arange(lo[1], hi[1] + 1) - 1.0) * dx
    a.pcolormesh(xs, ys, sub.T, cmap=cmap, norm=norm, edgecolors="0.6", linewidth=0.4, alpha=0.55)
    gi, gj = np.meshgrid(np.arange(lo[0], hi[0]), np.arange(lo[1], hi[1]), indexing="ij")
    flu = sub == FLUIDO
    a.plot((gi[flu] - 0.5) * dx, (gj[flu] - 0.5) * dx, "o", ms=3, color="#2a78d6", label="no de fluido")
    a.plot((gi[~flu] - 0.5) * dx, (gj[~flu] - 0.5) * dx, "s", ms=3, color="0.2", label="no solido")
    L = case.links
    if L.n.size:
        from lbm3d import CV
        pos = np.array(np.unravel_index(L.n, case.shape)).T
        c = CV[L.d]
        sel = (pos[:, ax] == idx) & (c[:, ax] == 0)
        sel &= (pos[:, h] >= lo[0]) & (pos[:, h] < hi[0]) & (pos[:, v] >= lo[1]) & (pos[:, v] < hi[1])
        if sel.any():
            q = q_links(case)[sel]
            p0 = (pos[sel][:, [h, v]] - 0.5)
            p1 = p0 + q[:, None] * c[sel][:, [h, v]]
            from matplotlib.collections import LineCollection
            a.add_collection(LineCollection(np.stack([p0, p1], 1) * dx, colors="#d03b3b", lw=0.8,
                                            label="link ate a parede (q)"))
            a.plot(p1[:, 0] * dx, p1[:, 1] * dx, "x", ms=3.5, color="#d03b3b", label="ponto de parede")
    a.set_aspect("equal")
    a.set_xlim(xs[0], xs[-1])
    a.set_ylim(ys[0], ys[-1])
    a.set_xlabel(_rotulo(case, h))
    a.set_ylabel(_rotulo(case, v))
    a.set_title("Detalhe da malha junto a parede (corte %s = no %d)" % (EIXOS[ax], idx), fontsize=9)
    a.legend(loc="upper center", bbox_to_anchor=(0.5, -0.09), ncol=2, fontsize=7, frameon=False)
    arquivos.append(_salvar(fig, os.path.join(pasta, "malha_detalhe.png")))

    # 3) resumo em texto
    n = int(np.prod(case.shape))
    nf = int((~case.solid).sum())
    q = q_links(case)
    curva = L.tag == 1
    linhas = ["Malha - caso %s" % case.name, "=" * 60,
              "Rede de velocidades     : D3Q27 (27 velocidades discretas, c_s^2 = 1/3)",
              "Malha cartesiana        : uniforme, dx = dt = 1 (unidades de rede)",
              "Dimensoes (nx, ny, nz)  : %d x %d x %d" % tuple(case.shape),
              "Nos totais              : %d" % n,
              "Nos de fluido           : %d (%.1f %%)" % (nf, 100.0 * nf / n),
              "Nos solidos             : %d" % (n - nf),
              "Links fluido->solido    : %d" % L.n.size]
    if curva.any():
        linhas.append("  em paredes curvas    : %d  (q min / med / max = %.3f / %.3f / %.3f)"
                      % (curva.sum(), q[curva].min(), q[curva].mean(), q[curva].max()))
    if getattr(case, "dx_fisico", None):
        linhas.append("Espacamento fisico      : dx = %.6g %s" % (case.dx_fisico, case.unidade))
    linhas.append("Viscosidade de rede     : nu = %.6g  (omega = %.5f, tau = %.5f)"
                  % (case.nu, case.omega, 1.0 / case.omega))
    if prm is not None:
        linhas.append("Modelo de colisao       : %s" % prm["model"])
    linhas += [""] + list(getattr(case, "descricao_malha", lambda: [])())
    caminho = os.path.join(pasta, "malha.txt")
    with open(caminho, "w", encoding="utf-8") as fh:
        fh.write("\n".join(linhas) + "\n")
    arquivos.append(caminho)
    return arquivos


def preprocessamento(case, pastas, prm=None):
    """Gera as imagens do dominio, das condicoes de contorno e da malha."""
    arquivos = figuras_dominio(case, pastas.dominio)
    arquivos += figuras_contorno(case, pastas.contorno)
    arquivos += figuras_malha(case, pastas.malha, prm)
    return arquivos


# ---------------------------------------------------------------------------
# Escoamento
# ---------------------------------------------------------------------------
def figuras_escoamento(sim, pasta):
    case = sim.case
    fluido = sim.fluid
    u = np.where(fluido[None], sim.u, np.nan)
    umag = np.sqrt((u ** 2).sum(axis=0))
    p = np.where(fluido, (sim.rho - 1.0) / 3.0, np.nan)  # p - p0 = c_s^2 (rho - 1)
    arquivos = []
    titulo = "%s - iteracao %d" % (case.name, sim.it)

    # |u| nos tres cortes
    fig, eixos = _grade_cortes(case)
    vmax = np.nanmax(umag) if np.isfinite(umag).any() else 1.0
    for a, (ax, idx) in zip(eixos, planos(case)):
        B, h, v = _corte(umag, ax, idx)
        im = a.imshow(B.T, origin="lower", extent=_extent(case, h, v), cmap="viridis",
                      vmin=0, vmax=vmax, interpolation="nearest", aspect="equal")
        a.set_xlabel(_rotulo(case, h))
        a.set_ylabel(_rotulo(case, v))
        a.set_title("|u|, corte %s = no %d" % (EIXOS[ax], idx), fontsize=9)
        fig.colorbar(im, ax=a, shrink=0.85, label="|u| [rede]")
    fig.suptitle(titulo)
    arquivos.append(_salvar(fig, os.path.join(pasta, "velocidade_cortes.png")))

    # componentes e pressao no corte principal
    ax, idx = planos(case)[0]
    fig = _figura(14, 8)
    campos = [("u_x", u[0], "RdBu_r"), ("u_y", u[1], "RdBu_r"), ("u_z", u[2], "RdBu_r"),
              ("p - p0 = (rho - 1)/3", p, "coolwarm")]
    for n, (nome, F, cm) in enumerate(campos):
        a = fig.add_subplot(2, 2, n + 1)
        B, h, v = _corte(F, ax, idx)
        lim = np.nanmax(np.abs(B)) if np.isfinite(B).any() else 1.0
        lim = lim if lim > 0 else 1.0
        im = a.imshow(B.T, origin="lower", extent=_extent(case, h, v), cmap=cm, vmin=-lim, vmax=lim,
                      interpolation="nearest", aspect="equal")
        fig.colorbar(im, ax=a, shrink=0.85)
        a.set_xlabel(_rotulo(case, h))
        a.set_ylabel(_rotulo(case, v))
        a.set_title(nome, fontsize=9)
    fig.subplots_adjust(hspace=0.35)
    fig.suptitle("%s  |  corte %s = no %d" % (titulo, EIXOS[ax], idx))
    arquivos.append(_salvar(fig, os.path.join(pasta, "componentes_pressao.png")))

    # linhas de corrente no corte principal (se houver escoamento no plano)
    B, h, v = _corte(umag, ax, idx)
    Uh, _, _ = _corte(np.nan_to_num(u[h]), ax, idx)
    Uv, _, _ = _corte(np.nan_to_num(u[v]), ax, idx)
    if np.nanmax(np.hypot(Uh, Uv)) > 1e-3 * max(vmax, 1e-300) and min(B.shape) >= 3:
        dx = getattr(case, "dx_fisico", None) or 1.0
        xs = (np.arange(B.shape[0]) - 0.5) * dx
        ys = (np.arange(B.shape[1]) - 0.5) * dx
        fig = _figura(12, 5.5)
        a = fig.add_subplot(1, 1, 1)
        im = a.imshow(B.T, origin="lower", extent=_extent(case, h, v), cmap="viridis",
                      interpolation="nearest", aspect="equal")
        try:
            a.streamplot(xs, ys, Uh.T, Uv.T, color="w", linewidth=0.6, density=1.4, arrowsize=0.6)
        except Exception:  # pragma: no cover - streamplot pode falhar em campos degenerados
            pass
        a.set_xlim(xs[0] - 0.5 * dx, xs[-1] + 0.5 * dx)
        a.set_ylim(ys[0] - 0.5 * dx, ys[-1] + 0.5 * dx)
        fig.colorbar(im, ax=a, shrink=0.8, label="|u| [rede]")
        a.set_xlabel(_rotulo(case, h))
        a.set_ylabel(_rotulo(case, v))
        a.set_title("Linhas de corrente (componentes no plano) - %s" % titulo, fontsize=9)
        arquivos.append(_salvar(fig, os.path.join(pasta, "linhas_corrente.png")))
    return arquivos


# ---------------------------------------------------------------------------
# Graficos de convergencia
# ---------------------------------------------------------------------------
def _serie(hist, k):
    """(iteracoes, valores) alinhados pelo fim; listas vazias ou de tamanhos diferentes
    nunca geram IndexError."""
    it = np.asarray(hist.get("it", []), float)
    v = np.asarray(hist.get(k, []), float)
    n = min(it.size, v.size)
    return it[it.size - n:], v[v.size - n:]


def graficos_convergencia(sim, pasta, tol=None):
    h = sim.hist
    arquivos = []
    if not h.get("it") or not h.get("res"):
        return arquivos
    fig = _figura(9, 5)
    a = fig.add_subplot(1, 1, 1)
    estilos = {"res": ("k", 2.0, "global ||du|| / ||u||"), "res_ux": ("#2a78d6", 1.0, "u_x"),
               "res_uy": ("#1baf7a", 1.0, "u_y"), "res_uz": ("#b87900", 1.0, "u_z"),
               "res_rho": ("#d03b3b", 1.0, "rho")}
    for k, (cor, lw, rot) in estilos.items():
        if k in h and len(h[k]):
            it, v = _serie(h, k)
            ok = v > 0
            if ok.any():
                a.semilogy(it[ok], v[ok], color=cor, lw=lw, label=rot)
    if tol and sim.case.steady:
        a.axhline(tol, color="r", ls="--", lw=1, label="tolerancia %.0e" % tol)
    a.grid(True, which="both", alpha=0.3)
    a.set_xlabel("iteracao")
    a.set_ylabel("residuo (variacao em %d it.)" % getattr(sim.args, "report_every", 100))
    a.legend(fontsize=8)
    a.set_title("Residuos - %s (it = %d, res = %.3e)" % (sim.case.name, sim.it, h["res"][-1]))
    arquivos.append(_salvar(fig, os.path.join(pasta, "convergencia.png")))

    extras = [k for k in h if k != "it" and not k.startswith("res") and _serie(h, k)[1].size]
    if extras:
        refs = sim.case.ref_ranges()
        fig = _figura(9, 2.4 * len(extras) + 0.8)
        for n, k in enumerate(extras):
            a = fig.add_subplot(len(extras), 1, n + 1)
            it, v = _serie(h, k)
            a.plot(it, v, "k.-", ms=3, lw=0.8)
            if k in refs:
                lo, hi = refs[k]
                if hi > lo:
                    a.axhspan(lo, hi, color="#1baf7a", alpha=0.18, label="faixa de referencia")
                else:
                    a.axhline(lo, color="#1baf7a", ls="--", label="referencia")
                a.legend(fontsize=7, loc="best")
            a.set_ylabel(k)
            a.grid(True, alpha=0.3)
            a.set_title("%s = %.6g" % (k, v[-1]), fontsize=8)
        fig.axes[-1].set_xlabel("iteracao")
        fig.suptitle("Grandezas monitoradas - %s" % sim.case.name)
        arquivos.append(_salvar(fig, os.path.join(pasta, "grandezas.png")))
    return arquivos


# ---------------------------------------------------------------------------
# Transferencia de calor (caso calor, [CN26])
# ---------------------------------------------------------------------------
def _graus(theta):
    """theta = (T - Tw)/(Tinf - Tw) -> graus C de [CN26] (Tw = 15 C, Tinf = 35 C)."""
    return 15.0 + 20.0 * theta


def figuras_temperatura(sim, pasta):
    """Campo de temperatura em graus C, janela -4 <= x/D <= 12, -4 <= y/D <= 4 (Figs. 4, 8, 11)."""
    case = sim.case
    th = sim.termico
    T = _graus(th.temperatura()[:, :, 0])
    if case.metodo == "sem_solido":
        T = np.where(case.corpo[:, :, 0], np.nan, T)
    D = case.D
    x = (np.arange(T.shape[0]) - 0.5 - case.xc) / D
    y = (np.arange(T.shape[1]) - 0.5 - case.yc) / D
    sx = (x >= -4) & (x <= 12)
    sy = (y >= -4) & (y <= 4)
    fig = _figura(9, 4.4)
    a = fig.add_subplot(1, 1, 1)
    import matplotlib
    cmap = matplotlib.colormaps["jet"].copy()
    cmap.set_bad("white")
    im = a.imshow(T[np.ix_(sx, sy)].T, origin="lower", cmap=cmap, vmin=15, vmax=35, aspect="equal",
                  extent=(x[sx][0] - 0.5 / D, x[sx][-1] + 0.5 / D, y[sy][0] - 0.5 / D, y[sy][-1] + 0.5 / D),
                  interpolation="bilinear")
    fig.colorbar(im, ax=a, shrink=0.8, label="T [°C]")
    a.set_xlabel("x/D")
    a.set_ylabel("y/D")
    from lbm3d import NOME_FORMA
    a.set_title("%s - Re = %g - %s  (t/Δt = %d)" % (NOME_FORMA[case.forma], case.Re,
                                                     case.metodo.replace("_", " "), sim.it), fontsize=10)
    arqs = [_salvar(fig, os.path.join(pasta, "temperatura.png"), dpi=130)]
    return arqs


def graficos_calor(case, sim, pasta, r):
    """Graficos de Nu no formato de [CN26]."""
    from lbm3d import NOME_FORMA, ciclos_completos, nu_churchill_bernstein
    arqs = []
    nome = "%s, Re = %g" % (NOME_FORMA[case.forma], case.Re)
    if case.metodo == "sem_solido":
        nu = np.asarray(case.serie_nu, float)
        n = nu.size
        if n < 2:
            return arqs
        inst = np.arange(1, n + 1)
        Nu = r["Nu"]
        # (a) historia completa + media
        fig = _figura(9, 3.6)
        a = fig.add_subplot(1, 1, 1)
        a.plot(inst, nu, "k-", lw=0.6, label="Nu_inst")
        a.axhline(Nu, color="b", lw=1.5, label="Nu = %.4f" % Nu)
        if case.forma == "circulo":
            a.axhline(nu_churchill_bernstein(case.Re, case.Pr), color="r", ls="--", lw=1,
                      label="Eq. 17 = %.3f" % nu_churchill_bernstein(case.Re, case.Pr))
        lo, hi = np.percentile(nu[n // 2:], [0, 100])
        pad = max(0.15 * (hi - lo), 0.05 * abs(Nu))
        a.set_ylim(lo - pad, hi + pad)
        a.set_xlabel("inst (t/Δt)")
        a.set_ylabel("Nu_inst")
        a.grid(alpha=0.3)
        a.legend(fontsize=8)
        a.set_title("Variação do número de Nusselt - %s (sem sólido)" % nome, fontsize=10)
        arqs.append(_salvar(fig, os.path.join(pasta, "nusselt_historia.png")))
        if case.transiente:
            ini, per = ciclos_completos(nu, 9)
            if ini is not None and per:
                # Fig. 10: sete periodos finais
                m = int(round(7 * per))
                fig = _figura(7, 3.4)
                a = fig.add_subplot(1, 1, 1)
                a.plot(np.arange(1, m + 1), nu[-m:], "k-", lw=1, label="Nu_inst")
                a.axhline(Nu, color="b", lw=2, label="Nu=%.4f" % Nu)
                a.set_xlabel("inst (t/Δt)")
                a.set_ylabel("Nu_inst")
                a.grid(alpha=0.3)
                a.legend(fontsize=8, loc="upper center", ncol=2)
                a.set_title(NOME_FORMA[case.forma], fontsize=11, fontweight="bold")
                arqs.append(_salvar(fig, os.path.join(pasta, "nusselt_7periodos.png")))
                # Fig. 5(b, c): sinal Nu_inst - Nu_avg nos 9 ciclos finais e FFT
                sinal = nu[ini:] - nu[ini:].mean()
                F = np.abs(np.fft.rfft(sinal))
                tempo = np.arange(sinal.size)       # em passos de rede
                fig = _figura(11, 3.6)
                a = fig.add_subplot(1, 2, 1)
                a.plot(tempo * case.U / case.D, sinal, "k-", lw=0.8)
                a.set_xlabel("t U/D (amostra)")
                a.set_ylabel("Nu_inst - Nu_avg")
                a.grid(alpha=0.3)
                a.set_title("(b) amostra: 9 ciclos completos", fontsize=9)
                a = fig.add_subplot(1, 2, 2)
                a.plot(np.arange(F.size)[:40], F[:40], "b.-")
                a.set_xlabel("frequência (ciclos na amostra)")
                a.set_ylabel("Magnitude |FFT|")
                a.grid(alpha=0.3)
                a.set_title("(c) pico em %d ciclos; St = D/(T U) = %.4f"
                            % (int(np.argmax(F)), case.D / (per * case.U)), fontsize=9)
                arqs.append(_salvar(fig, os.path.join(pasta, "nusselt_fft.png")))
        return arqs
    # --- metodo com solido: Figs. 6 e 7 ---
    tc = np.asarray(case.serie_tc, float)
    n = tc.size
    if n < 20:
        return arqs
    t = np.arange(n, dtype=float)
    fator = case.perimetro / (case.Cs * case.area_discreta)
    h_ref = nu_churchill_bernstein(case.Re, case.Pr) * case.alfa / case.D if case.forma == "circulo" else None
    h = r.get("h_rede", np.nan)
    fig = _figura(12, 4)
    a = fig.add_subplot(1, 2, 1)
    a.plot(t, _graus(tc), "k-", lw=2, label="CFD (centro do sólido)")
    if h_ref:
        a.plot(t, _graus(1 - np.exp(-h_ref * fator * t)), "r--", lw=1.2, label="Eq. 14 com h da Eq. 17")
    a.plot(t, _graus(1 - np.exp(-h * fator * t)), "b:", lw=1.2, label="Eq. 14 com h ajustado")
    a.set_xlabel("inst (t/Δt)")
    a.set_ylabel("T [°C]")
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
    a.set_title("Perfil de temperatura do sólido - %s" % nome, fontsize=9)
    a = fig.add_subplot(1, 2, 2)
    pontos = np.unique(np.linspace(max(20, n // 50), n, 40).astype(int))
    nus = []
    for m in pontos:
        c0 = case.serie_tc
        case.serie_tc = c0[:m]
        nus.append(case.ajuste_h(0.0)[0] * case.D / case.alfa)
        case.serie_tc = c0
    a.plot(pontos, nus, "k.-")
    if case.forma == "circulo":
        a.axhline(nu_churchill_bernstein(case.Re, case.Pr), color="r", ls="--", label="Eq. 17")
        a.legend(fontsize=8)
    a.set_xlabel("número acumulado de passos")
    a.set_ylabel("Nu (ajuste com todos os passos)")
    a.grid(alpha=0.3)
    a.set_title("Nu em função dos passos acumulados", fontsize=9)
    arqs.append(_salvar(fig, os.path.join(pasta, "solido_perfil_nusselt.png")))
    # Fig. 7: soma dos quadrados x h
    _, (ii, tt, Tcfd, soma) = case.ajuste_h()
    hs = np.linspace(0.8 * h, 1.2 * h, 81)
    fig = _figura(6.5, 4)
    a = fig.add_subplot(1, 1, 1)
    a.plot(hs * case.D / case.alfa, [soma(v) for v in hs], "k-")
    a.axvline(h * case.D / case.alfa, color="b", ls="--", label="mínimo: Nu = %.4f" % (h * case.D / case.alfa))
    a.set_xlabel("Nu = h D / k_f")
    a.set_ylabel("Σ dif²  (inst = %d ... %d)" % (ii + 1, n))
    a.grid(alpha=0.3)
    a.legend(fontsize=8)
    a.set_title("Determinação do coeficiente h - %s" % nome, fontsize=9)
    arqs.append(_salvar(fig, os.path.join(pasta, "solido_ajuste_h.png")))
    return arqs
