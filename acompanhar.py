#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
acompanhar.py - acompanha a convergencia de uma simulacao do lbm3d.py em andamento
(ou terminada) lendo <pasta>/convergencia.csv, gravado a cada 100 iteracoes.

Funciona com simulacoes rodando sem janela (--no-gui), em segundo plano ou em outra
maquina (basta ter acesso a pasta).  Mostra:
  * curva de residuo ||u^n - u^(n-100)|| / ||u^n|| com a tolerancia de parada;
  * taxa de convergencia (decadas por 1000 it.), velocidade (it/s) e estimativa
    do tempo restante ate atingir a tolerancia (extrapolacao log-linear);
  * evolucao das grandezas do caso (C_D, C_L, C, erro L2 ...), com a variacao
    relativa nos ultimos relatorios.

Uso:
  python acompanhar.py saida_spheres                # gera e abre o painel HTML
  python acompanhar.py saida_spheres --servir       # painel HTML ao vivo em http://127.0.0.1:8765
  python acompanhar.py saida_cylinder --janela      # janela matplotlib atualizada a cada 5 s
  python acompanhar.py saida_cylinder --texto       # so terminal (ex.: via ssh)
  python acompanhar.py saida_cylinder --png         # grava acompanhamento.png e sai

O lbm3d.py ja grava <saida>/convergencia.html a cada 100 iteracoes: basta abrir esse
arquivo no navegador (ele se recarrega sozinho).
"""
import argparse
import os
import sys
import time

import numpy as np


def ler(pasta):
    path = os.path.join(pasta, "convergencia.csv")
    if not os.path.exists(path):
        return None, None
    meta = {}
    with open(path) as fh:
        first = fh.readline()
        if first.startswith("#"):
            for tok in first[1:].split():
                if "=" in tok:
                    k, v = tok.split("=", 1)
                    meta[k] = v
            header = fh.readline().strip().split(",")
        else:
            header = first.strip().split(",")
        linhas = [l.strip().split(",") for l in fh if l.strip() and not l.startswith(("#", "it,"))]
    if not linhas:
        return meta, None
    dados = np.array([[float(x) for x in l] for l in linhas if len(l) == len(header)])
    cols = {h: dados[:, i] for i, h in enumerate(header)}
    # se houve reinicio a partir de um checkpoint anterior, fica so a ultima passagem
    it = cols["it"]
    quebra = np.where(np.diff(it) <= 0)[0]
    if quebra.size:
        ini = quebra[-1] + 1
        cols = {k: v[ini:] for k, v in cols.items()}
    return meta, cols


def diagnostico(meta, c, janela=20):
    it, res, t = c["it"], c["res"], c["tempo_unix"]
    tol = float(meta.get("tol", 1e-7))
    out = {"it": int(it[-1]), "res": res[-1], "tol": tol, "taxa": np.nan,
           "its": np.nan, "eta_it": np.nan, "eta_s": np.nan}
    n = min(janela, len(it))
    if n >= 3:
        ok = res[-n:] > 0
        if ok.sum() >= 3:
            slope = np.polyfit(it[-n:][ok], np.log10(res[-n:][ok]), 1)[0]
            out["taxa"] = -slope * 1000.0  # decadas por 1000 iteracoes
            dt = t[-1] - t[-n]
            if dt > 0:
                out["its"] = (it[-1] - it[-n]) / dt
            if slope < 0 and res[-1] > tol:
                out["eta_it"] = (np.log10(tol) - np.log10(res[-1])) / slope
                out["eta_s"] = out["eta_it"] / out["its"] if out["its"] > 0 else np.nan
    out["ativo"] = (time.time() - t[-1]) < max(120.0, 5 * 100 / out["its"]) if out["its"] > 0 else False
    return out


def fmt_tempo(s):
    if not np.isfinite(s):
        return "?"
    s = int(s)
    return "%dh%02dm%02ds" % (s // 3600, (s % 3600) // 60, s % 60)


def texto_status(meta, c, d):
    extras = [k for k in c if k not in ("it", "res", "tempo_unix", "mlups")]
    partes = ["[%s] it=%d  res=%.3e (tol %.0e)" % (meta.get("caso", "?"), d["it"], d["res"], d["tol"])]
    if d["res"] <= d["tol"]:
        partes.append("CONVERGIDO")
    else:
        partes.append("taxa=%.3f dec/1000it  %.2f it/s  faltam ~%s it (%s)"
                      % (d["taxa"], d["its"], "?" if not np.isfinite(d["eta_it"]) else "%d" % d["eta_it"],
                         fmt_tempo(d["eta_s"])))
    for k in extras:
        v = c[k]
        rel = abs(v[-1] - v[-6]) / max(abs(v[-1]), 1e-300) if len(v) > 5 else np.nan
        partes.append("%s=%.6g (var. 500it: %.1e)" % (k, v[-1], rel))
    if not d["ativo"]:
        partes.append("(sem atualizacao recente: simulacao parada ou pausada)")
    return "  |  ".join(partes)


class Painel:
    def __init__(self, png=False):
        import matplotlib
        if png:
            matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        self.plt = plt
        self.fig = None

    def desenhar(self, meta, c, d, pasta):
        plt = self.plt
        extras = [k for k in c if k not in ("it", "res", "tempo_unix", "mlups")]
        n = 1 + max(1, len(extras))
        if self.fig is None or len(self.fig.axes) != n:
            if self.fig is not None:
                plt.close(self.fig)
            self.fig, axs = plt.subplots(n, 1, figsize=(9, 2.6 + 2.0 * n), sharex=True, squeeze=False)
        axs = self.fig.axes
        for ax in axs:
            ax.clear()
            ax.grid(True, alpha=0.3, which="both")
        it = c["it"]
        ax = axs[0]
        ax.semilogy(it, c["res"], "b.-", ms=3, label="residuo")
        ax.axhline(d["tol"], color="r", ls="--", label="tolerancia %.0e" % d["tol"])
        if np.isfinite(d["eta_it"]) and d["res"] > d["tol"]:
            ax.semilogy([it[-1], it[-1] + d["eta_it"]], [d["res"], d["tol"]], "b:", alpha=0.6,
                        label="extrapolacao")
        ax.set_ylabel(r"$\|u^n-u^{n-100}\|/\|u^n\|$")
        ax.legend(loc="upper right", fontsize=8)
        estado = ("CONVERGIDO" if d["res"] <= d["tol"] else
                  "faltam ~%s (%s)" % ("?" if not np.isfinite(d["eta_it"]) else "%d it" % d["eta_it"],
                                       fmt_tempo(d["eta_s"])))
        ax.set_title("%s  [%s]  it=%d  res=%.2e  |  %.2f it/s  |  %s%s"
                     % (pasta, meta.get("caso", "?"), d["it"], d["res"], d["its"], estado,
                        "" if d["ativo"] else "  |  PARADA/PAUSADA"), fontsize=9)
        for ax, k in zip(axs[1:], extras):
            ax.plot(it, c[k], "k.-", ms=3)
            ax.set_ylabel(k)
            v = c[k][-max(5, len(it) // 3):]
            lo, hi = np.min(v), np.max(v)
            pad = 0.1 * (hi - lo) + 1e-12 * max(abs(hi), 1.0)
            ax.set_ylim(lo - pad, hi + pad)
            ax.set_title("%s = %.6g  (zoom no ultimo terco)" % (k, c[k][-1]), fontsize=8)
        axs[-1].set_xlabel("iteracao")
        self.fig.tight_layout()


TEMPLATE = os.path.join(os.path.dirname(os.path.abspath(__file__)), "painel_template.html")


def gerar_html(pasta, intervalo=5):
    """Grava <pasta>/convergencia.html: pagina autocontida (sem internet) com os dados
    de convergencia.csv embutidos.  Aberta do disco ela se recarrega a cada `intervalo`
    segundos; servida por HTTP (--servir) ela busca o CSV sem recarregar."""
    import json
    csv = os.path.join(pasta, "convergencia.csv")
    texto = open(csv).read() if os.path.exists(csv) else ""
    with open(TEMPLATE, encoding="utf-8") as fh:
        tpl = fh.read()
    ini, fim = tpl.index("/*__DADOS__*/"), tpl.index("/*__FIM__*/")
    dados = json.dumps({"csv": texto, "intervalo": intervalo,
                        "pasta": os.path.basename(os.path.abspath(pasta))}).replace("</", "<\\/")
    html = tpl[:ini] + "/*__DADOS__*/" + dados + tpl[fim:]
    out = os.path.join(pasta, "convergencia.html")
    tmp = out + ".tmp"
    with open(tmp, "w", encoding="utf-8") as fh:
        fh.write(html)
    os.replace(tmp, out)
    return out


def servir(pasta, porta, intervalo):
    """Servidor HTTP local (127.0.0.1) para a pasta: a pagina busca o CSV ao vivo."""
    import functools
    import http.server
    import webbrowser

    class Quieto(http.server.SimpleHTTPRequestHandler):
        def log_message(self, *a):
            pass

        def end_headers(self):
            self.send_header("Cache-Control", "no-store")
            super().end_headers()

    gerar_html(pasta, intervalo)
    handler = functools.partial(Quieto, directory=os.path.abspath(pasta))
    srv = http.server.ThreadingHTTPServer(("127.0.0.1", porta), handler)
    url = "http://127.0.0.1:%d/convergencia.html" % srv.server_address[1]
    print("Painel ao vivo em %s  (Ctrl+C para encerrar)" % url)
    webbrowser.open(url)
    try:
        srv.serve_forever()
    except KeyboardInterrupt:
        pass


def main():
    p = argparse.ArgumentParser(description="Acompanha a convergencia de uma simulacao lbm3d.py")
    p.add_argument("pasta", help="pasta de saida da simulacao (contem convergencia.csv)")
    p.add_argument("--intervalo", type=float, default=5.0, help="segundos entre atualizacoes")
    p.add_argument("--html", action="store_true",
                   help="gera convergencia.html e abre no navegador (padrao)")
    p.add_argument("--servir", nargs="?", type=int, const=8765, default=None, metavar="PORTA",
                   help="serve a pasta em http://127.0.0.1:PORTA com dados ao vivo")
    p.add_argument("--janela", action="store_true", help="janela matplotlib em vez de HTML")
    p.add_argument("--texto", action="store_true", help="so imprime o estado no terminal")
    p.add_argument("--png", action="store_true", help="grava acompanhamento.png e sai")
    a = p.parse_args()

    if a.servir is not None:
        return servir(a.pasta, a.servir, a.intervalo)
    if not (a.janela or a.texto or a.png):
        import webbrowser
        out = gerar_html(a.pasta, a.intervalo)
        print("Painel gerado:", out)
        webbrowser.open("file://" + os.path.abspath(out))
        return

    painel = None if a.texto else Painel(png=a.png)
    ultimo = None
    try:
        while True:
            meta, c = ler(a.pasta)
            if c is None or len(c["it"]) < 1:
                print("aguardando %s/convergencia.csv ..." % a.pasta)
            elif ultimo != c["it"][-1]:
                ultimo = c["it"][-1]
                d = diagnostico(meta, c)
                print(time.strftime("%H:%M:%S"), texto_status(meta, c, d), flush=True)
                if painel is not None:
                    painel.desenhar(meta, c, d, a.pasta)
                    if a.png:
                        out = os.path.join(a.pasta, "acompanhamento.png")
                        painel.fig.savefig(out, dpi=110)
                        print("grafico salvo em", out)
                        return
            if a.png:
                return
            if painel is not None:
                painel.plt.pause(a.intervalo)
            else:
                time.sleep(a.intervalo)
    except KeyboardInterrupt:
        pass


if __name__ == "__main__":
    sys.exit(main())
