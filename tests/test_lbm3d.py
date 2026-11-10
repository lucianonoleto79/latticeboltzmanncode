# -*- coding: utf-8 -*-
"""Testes rapidos (poucos minutos) rodados a cada push pelo GitHub Actions.

Os benchmarks completos (cilindro de Schaefer-Turek, arranjo de esferas) demoram
horas e rodam pelo fluxo "Simulacao LBM" (.github/workflows/simulacao.yml).
"""
import json
import os
import re
import threading
import time

import numpy as np
import pytest

import acompanhar
import lbm3d


def rodar(tmp_path, *args):
    out = str(tmp_path)
    sim = lbm3d.main(["--no-gui", "--save-every", "0", "--no-figs", "--outdir", out, *args])
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
    # --restart aceita a pasta do caso (o checkpoint fica em checkpoint/checkpoint.npz)
    cont = lbm3d.main(["--no-gui", "--restart", out, "--max-iter", "400", "--save-every", "0",
                       "--no-figs"])
    assert cont.it == 400
    assert np.array_equal(direto.f, cont.f)


def test_limite_de_tempo(tmp_path):
    """--max-time encerra salvando o checkpoint, sem rodar a validacao final."""
    _, out = rodar(tmp_path, "--case", "duct", "--res", "12", "--max-time", "1e-6")
    estado = json.load(open(os.path.join(out, "estado.json")))
    assert not estado["terminou"] and estado["it"] == 100
    assert os.path.exists(os.path.join(out, "checkpoint", "checkpoint.npz"))
    assert not os.path.exists(os.path.join(out, "validacao.txt"))


def test_painel_html(tmp_path):
    _, out = rodar(tmp_path, "--case", "duct", "--res", "12", "--max-iter", "300")
    html = open(acompanhar.gerar_html(out), encoding="utf-8").read()
    assert "L2_error" in html and "/*__FIM__*/" in html
    meta, c = acompanhar.ler(out)
    assert meta["caso"] == "duct" and list(c["it"]) == [100, 200, 300]
    assert os.path.exists(os.path.join(out, "graficos", "convergencia.html"))
    assert {"res_ux", "res_uy", "res_uz", "res_rho"} <= set(c)


def test_pastas_do_caso(tmp_path):
    """Cada caso tem dominio/, contorno/, malha/, escoamento/, graficos/ e checkpoint/."""
    out = str(tmp_path / "caso")
    sim = lbm3d.main(["--no-gui", "--case", "duct", "--res", "12", "--max-iter", "200",
                      "--save-every", "100", "--outdir", out])
    esperado = ["dominio/dominio_3d.png", "dominio/dominio_cortes.png",
                "contorno/contorno_cortes.png", "contorno/contorno_3d.png",
                "contorno/condicoes_contorno.txt", "malha/malha_corte.png",
                "malha/malha_detalhe.png", "malha/malha.txt",
                "escoamento/parcial_it0000100/velocidade_cortes.png",
                "escoamento/parcial_it0000100/duct_it0000100.vtk",
                "escoamento/final/componentes_pressao.png", "escoamento/final/duct_it0000200.vtk",
                "graficos/convergencia.png", "graficos/grandezas.png", "graficos/duct_validacao.png",
                "graficos/convergencia.csv", "checkpoint/checkpoint.npz",
                "parametros.json", "execucao.log", "estado.json", "validacao.txt", "LEIAME.txt"]
    faltando = [f for f in esperado if not os.path.exists(os.path.join(out, f))]
    assert not faltando, faltando
    # residuos por componente: res^2 = res_ux^2 + res_uy^2 + res_uz^2
    r = sim.residuos()
    assert np.isclose(r["res"] ** 2, r["res_ux"] ** 2 + r["res_uy"] ** 2 + r["res_uz"] ** 2)
    log = open(os.path.join(out, "execucao.log"), encoding="utf-8").read()
    assert "res(ux)" in log and "res(rho)" in log


def test_mapa_contorno_cilindro():
    import figuras
    args = lbm3d.parse_args(["--case", "cylinder", "--res", "4"])
    case = lbm3d.Cylinder(args)
    cod = figuras.mapa_contorno(case)
    assert (cod[0, 2:-2, 2:-2] == figuras.MOVEL).all()          # entrada
    assert (cod[-1, 1:-1, 1:-1] == figuras.SAIDA).all()          # saida
    assert (cod == figuras.CURVA).any() and (cod == figuras.PAREDE).any()
    assert case.tag() == "cylinder_D4_Re20"


def test_arquivos_de_controle(tmp_path):
    """SALVAR grava um resultado parcial sem parar; PARAR salva o checkpoint e encerra."""
    out = str(tmp_path)
    res = {}
    t = threading.Thread(target=lambda: res.update(sim=lbm3d.main(
        ["--no-gui", "--case", "duct", "--res", "12", "--max-iter", "1000000", "--save-every", "0",
         "--no-figs", "--outdir", out])))
    t.start()
    try:
        csv = os.path.join(out, "graficos", "convergencia.csv")
        limite = time.time() + 60
        while not os.path.exists(csv) and time.time() < limite:
            time.sleep(0.05)
        assert acompanhar.controlar(out, "salvar") == 0
        while not any(d.startswith("parcial_") for d in os.listdir(os.path.join(out, "escoamento")))                 and time.time() < limite:
            time.sleep(0.05)
        assert acompanhar.controlar(out, "parar") == 0
        t.join(60)
    finally:
        if t.is_alive():  # garante o fim da thread se algo acima falhar
            open(os.path.join(out, "PARAR"), "w").close()
            t.join(60)
    assert not t.is_alive()
    estado = json.load(open(os.path.join(out, "estado.json")))
    assert not estado["terminou"] and estado["it"] < 1000000
    assert os.path.exists(os.path.join(out, "checkpoint", "checkpoint.npz"))
    assert not os.path.exists(os.path.join(out, "PARAR"))


def test_cilindro_transiente(tmp_path):
    """cylinder2z (3D-2Z, Re = 100): grava a serie temporal e analisa ao fim de --max-iter."""
    sim, out = rodar(tmp_path, "--case", "cylinder2z", "--res", "4", "--max-iter", "300")
    serie = np.loadtxt(os.path.join(out, "graficos", "serie_temporal.csv"), delimiter=",", skiprows=2)
    assert serie.shape == (60, 5) and serie[-1, 0] == 300
    txt = open(os.path.join(out, "validacao.txt")).read()
    assert "3D-2Z" in txt and "C_D max" in txt
    assert json.load(open(os.path.join(out, "estado.json")))["terminou"]
    assert sim.case.tag() == "cylinder2z_D4_Re100"
