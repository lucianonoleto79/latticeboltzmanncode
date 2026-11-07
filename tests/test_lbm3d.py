# -*- coding: utf-8 -*-
"""Testes rapidos (poucos minutos) rodados a cada push pelo GitHub Actions.

Os benchmarks completos (cilindro de Schaefer-Turek, arranjo de esferas) demoram
horas e rodam pelo fluxo "Simulacao LBM" (.github/workflows/simulacao.yml).
"""
import json
import os
import re

import numpy as np
import pytest

import acompanhar
import lbm3d


def rodar(tmp_path, *args):
    out = str(tmp_path)
    sim = lbm3d.main(["--no-gui", "--save-every", "0", "--outdir", out, *args])
    return sim, out


def erro_onda(out):
    txt = open(os.path.join(out, "validacao.txt")).read()
    return abs(float(re.search(r"erro relativo=([-+0-9.eE]+)", txt).group(1)))


def test_k17_simetrias():
    """Pesos somam 1 e as taxas K17 ficam em (0, 2) na faixa de uso (omega > 1.75;
    em omega = 1.75 a taxa w4 = 8(w-2)(4w-7)/(...) se anula)."""
    assert np.isclose(lbm3d.W.sum(), 1.0)
    assert lbm3d.collision_parameters(1.75, "k17", None)["model"] == "k17"
    for w in (1.8, 1.85, 1.95):
        prm = lbm3d.collision_parameters(w, "k17", None)
        assert prm["model"] == "k17"
        assert all(0.0 < prm[k] < 2.0 for k in ("w3", "w4", "w5"))
    # abaixo de 1.75 o codigo recai no cumulante classico
    assert lbm3d.collision_parameters(1.5, "k17", None)["model"] == "cumulant"


def test_onda_cisalhamento_k17_mais_precisa(tmp_path):
    """Viscosidade medida no decaimento da onda: K17 ~4a ordem, muito melhor que BGK [G17]."""
    args = ["--case", "shearwave", "--res", "16", "--nu", "0.01", "--u-adv", "0.05"]
    _, out_k = rodar(tmp_path / "k17", *args, "--model", "k17")
    _, out_b = rodar(tmp_path / "bgk", *args, "--model", "bgk")
    ek, eb = erro_onda(out_k), erro_onda(out_b)
    assert ek < 1e-3, ek
    assert eb > 10 * ek, (ek, eb)


def test_duto_retangular(tmp_path):
    """Perfil de Poiseuille em duto retangular contra a serie de White (2006), Eq. 3-48."""
    sim, out = rodar(tmp_path, "--case", "duct", "--res", "12", "--tol", "1e-7")
    erro = sim.case.error(sim)[0]
    assert erro < 1e-2, erro
    estado = json.load(open(os.path.join(out, "estado.json")))
    assert estado["convergiu"] and estado["terminou"]


def test_reinicio_bit_exato(tmp_path):
    """Pausar (checkpoint) e continuar reproduz exatamente a rodada sem interrupcao."""
    caso = ["--case", "duct", "--res", "12"]
    direto, _ = rodar(tmp_path / "direto", *caso, "--max-iter", "400")
    _, out = rodar(tmp_path / "partes", *caso, "--max-iter", "200")
    estado = json.load(open(os.path.join(out, "estado.json")))
    assert estado["it"] == 200
    cont = lbm3d.main(["--no-gui", "--restart", os.path.join(out, "checkpoint.npz"),
                       "--max-iter", "400", "--outdir", out, "--save-every", "0"])
    assert cont.it == 400
    assert np.array_equal(direto.f, cont.f)


def test_limite_de_tempo(tmp_path):
    """--max-time encerra salvando o checkpoint, sem rodar a validacao final."""
    _, out = rodar(tmp_path, "--case", "duct", "--res", "12", "--max-time", "1e-6")
    estado = json.load(open(os.path.join(out, "estado.json")))
    assert not estado["terminou"] and estado["it"] == 100
    assert os.path.exists(os.path.join(out, "checkpoint.npz"))
    assert not os.path.exists(os.path.join(out, "validacao.txt"))


def test_painel_html(tmp_path):
    _, out = rodar(tmp_path, "--case", "duct", "--res", "12", "--max-iter", "300")
    html = open(acompanhar.gerar_html(out), encoding="utf-8").read()
    assert "L2_error" in html and "/*__FIM__*/" in html
    meta, c = acompanhar.ler(out)
    assert meta["caso"] == "duct" and list(c["it"]) == [100, 200, 300]
