# -*- coding: utf-8 -*-
"""Testes do nucleo em C, da transferencia de calor (caso calor, [CN26]), do contador de
tempo de simulacao e da robustez contra 'list index out of range'."""
import json
import os

import numpy as np
import pytest

import acompanhar
import figuras
import lbm3d

C_DISPONIVEL = lbm3d.kernel_c() is not None


def rodar(tmp_path, *args):
    out = str(tmp_path)
    sim = lbm3d.main(["--no-gui", "--save-every", "0", "--no-figs", "--outdir", out, *args])
    return sim, out


# ---------------------------------------------------------------------------
# Nucleo em C
# ---------------------------------------------------------------------------
@pytest.mark.skipif(not C_DISPONIVEL, reason="sem compilador C")
@pytest.mark.parametrize("modelo,omega,lim", [("k17", 1.9, None), ("k17", 1.9, 0.01),
                                              ("cumulant", 1.3, None), ("bgk", 1.6, None)])
def test_nucleo_c_igual_numpy(modelo, omega, lim):
    """Colisao + propagacao em C = NumPy (arredondamento), com solidos e forca."""
    so = lbm3d.kernel_c()
    rng = np.random.default_rng(3)
    for shape in ((9, 7, 5), (6, 40, 1)):
        rho = 1 + 0.01 * rng.standard_normal(shape)
        u = [0.05 * rng.standard_normal(shape) for _ in range(3)]
        f = np.ascontiguousarray(lbm3d.equilibrium(rho, *u) * (1 + 0.02 * rng.standard_normal((3, 3, 3) + shape)))
        solid = rng.random(shape) < 0.15
        F = (1e-4, -2e-5, 3e-5)
        prm = lbm3d.collision_parameters(omega, modelo, lim)
        a = f.copy()
        lbm3d.collide(a, F, prm)
        a.reshape(27, -1)[:, np.flatnonzero(solid)] = lbm3d.W.reshape(27, 1)
        an = np.empty_like(a)
        lbm3d.stream(a, an)
        bn = np.empty_like(f)
        r, ux, uy, uz = (np.empty(shape) for _ in range(4))
        so.collide_stream(f, bn, solid.astype(np.uint8), *shape, np.array(F, float), lbm3d.prm_c(prm), r, ux, uy, uz)
        assert np.abs(an - bn).max() < 1e-13


@pytest.mark.skipif(not C_DISPONIVEL, reason="sem compilador C")
@pytest.mark.parametrize("argv", [["--case", "cylinder", "--res", "6"], ["--case", "spheres", "--res", "12"],
                                  ["--case", "calor", "--res", "6", "--Ds", "10", "--Re", "10"],
                                  ["--case", "calor", "--res", "6", "--Ds", "10", "--Re", "100",
                                   "--metodo", "com_solido"]])
def test_simulacao_c_igual_numpy(argv):
    """Passo completo (paredes, saida, forca, campo termico) igual com e sem o nucleo em C."""
    res = []
    for sem_c in (True, False):
        a = lbm3d.parse_args(argv + ["--no-gui"])
        a.sem_c = sem_c
        case = lbm3d.CASES[a.case](a)
        sim = lbm3d.Simulation(case, a)
        for _ in range(40):
            sim.step()
        res.append((sim.f.copy(), sim.body_force.copy(),
                    None if sim.termico is None else sim.termico.g.copy()))
    assert np.abs(res[0][0] - res[1][0]).max() < 1e-12
    assert np.allclose(res[0][1], res[1][1], rtol=1e-9, atol=1e-12)
    if res[0][2] is not None:
        assert np.abs(res[0][2] - res[1][2]).max() < 1e-11


# ---------------------------------------------------------------------------
# Transferencia de calor
# ---------------------------------------------------------------------------
class _Obj:
    pass


def test_conducao_anel_parede_curva():
    """Conducao entre cilindros concentricos (T = 0 em r1, T = 1 em r2): calor exato
    Q = 2 pi alfa / ln(r2/r1); paredes curvas pelo anti-bounce-back interpolado [LMK13]."""
    n, r1, r2, alfa = 60, 6.3, 25.0, 0.1
    sh = (n, n, 1)
    X, Y = np.meshgrid(np.arange(n) - 0.5, np.arange(n) - 0.5, indexing="ij")
    xc = yc = n / 2 - 0.5 + 0.23
    R = np.hypot(X - xc, Y - yc)[:, :, None]
    dentro, fora = R < r1, R > r2

    def qcirc(rad, externo):
        def q(x, y, z, cx, cy, cz):
            px, py = x - xc, y - yc
            a, b, c = cx * cx + cy * cy, px * cx + py * cy, px * px + py * py - rad ** 2
            d = np.sqrt(np.maximum(b * b - a * c, 0))
            return ((-b + d) if externo else (-b - d)) / a
        return q
    case = _Obj()
    case.shape, case.solid = sh, dentro | fora
    case.saida_termica = lambda g: None
    case.apos_passo_termico = lambda s, th: None
    sim = _Obj()
    sim.case, sim.u, sim.so = case, np.zeros((3,) + sh), lbm3d.kernel_c()
    th = lbm3d.Termico(sim, ~case.solid, np.ones(sh), np.ones(sh), np.full(sh, 0.5),
                       [dict(mask=dentro, Tw=0.0, q=qcirc(r1, False), tag=1),
                        dict(mask=fora, Tw=1.0, q=qcirc(r2, True), tag=2)], alfa)
    Q = []
    case.apos_passo_termico = lambda s, t: Q.append(t.links.calor(t.gtmp, 1))
    for _ in range(int(5 * r2 ** 2 / alfa)):
        th.passo(sim)
    exato = 2 * np.pi * alfa / np.log(r2 / r1)
    assert abs(Q[-1] - exato) / exato < 5e-3, (Q[-1], exato)


def test_nusselt_referencia_eq17():
    """Eq. 17 de [CN26] reproduz os valores de referencia do artigo (Tabela 4)."""
    for Re, ref in ((1, 0.791), (10, 1.853), (100, 5.210), (200, 7.244)):
        assert round(lbm3d.nu_churchill_bernstein(Re), 3) == ref


def test_formas_geometria():
    """Perimetros e areas das secoes da Fig. 9 de [CN26] (D = altura) conferem com a mascara."""
    D = 40
    x = np.arange(-1.5 * D, 1.5 * D, 0.25) + 0.125
    X, Y = np.meshgrid(x, x, indexing="ij")
    for forma in lbm3d.FORMAS:
        dentro, per, area = lbm3d.geometria_forma(forma, D)
        m = dentro(X, Y)
        assert abs(m.sum() * 0.25 ** 2 - area) / area < 0.02, forma
        alt = np.ptp(Y[m]) + 0.25
        assert abs(alt - D) / D < 0.02, forma
        assert per > 0


def test_calor_reinicio_bit_exato(tmp_path):
    """Reinicio do caso de calor (campo termico e series de Nu no checkpoint) e bit a bit."""
    base = ["--case", "calor", "--res", "6", "--Ds", "10", "--Re", "100", "--report-every", "10"]
    simA, _ = rodar(tmp_path / "a", *base, "--max-iter", "60")
    rodar(tmp_path / "b", *base, "--max-iter", "30")
    simB = lbm3d.main(["--no-gui", "--save-every", "0", "--no-figs", "--restart", str(tmp_path / "b"),
                       "--max-iter", "60"])
    assert np.array_equal(simA.f, simB.f)
    assert np.array_equal(simA.termico.g, simB.termico.g)
    assert simA.case.serie_nu == simB.case.serie_nu


def test_calor_criterio_estacionario_e_resultado(tmp_path):
    """Caso estacionario para pelo criterio da Eq. 16 e grava resultado_calor.json."""
    a = ["--case", "calor", "--res", "5", "--Ds", "8", "--Re", "10", "--report-every", "50"]
    a = lbm3d.parse_args(a + ["--no-gui"])
    case = lbm3d.CASES["calor"](a)
    case.tol_nu = 1e-3                       # criterio frouxo para o teste ser rapido
    lbm3d.CASES["_teste_calor"] = lambda args: case
    try:
        sim, out = rodar(tmp_path, "--case", "_teste_calor", "--report-every", "50", "--max-iter", "20000")
    finally:
        lbm3d.CASES.pop("_teste_calor")
    r = json.load(open(os.path.join(out, "resultado_calor.json")))
    assert "Eq. 16 adaptada" in r["parada"]
    assert 1.0 < r["Nu"] < 5.0                # dominio pequeno: so a ordem de grandeza
    est = json.load(open(os.path.join(out, "estado.json")))
    assert est["terminou"] and est["convergiu"]


def test_ciclos_completos_com_deriva():
    """Periodo do Nu_inst achado mesmo com deriva lenta (inicio do transiente)."""
    n = 16000
    t = np.arange(n)
    x = 5.0 + 0.2 * np.exp(-t / 3000.0) + 0.004 * np.sin(2 * np.pi * t / 410.0)
    ini, per = lbm3d.ciclos_completos(x, 9)
    assert ini is not None and abs(per - 410) < 5
    assert abs((n - ini) - 9 * per) <= 1
    assert lbm3d.ciclos_completos(np.full(n, 3.0))[0] is None


def test_com_solido_ajuste_h():
    """O ajuste de h por minimos quadrados recupera h de uma curva exponencial exata."""
    a = lbm3d.parse_args(["--case", "calor", "--res", "6", "--Ds", "8", "--Re", "10",
                          "--metodo", "com_solido", "--no-gui"])
    case = lbm3d.CASES["calor"](a)
    h = 0.0123
    fator = case.perimetro / (case.Cs * case.area_discreta)
    case.serie_tc = list(1 - np.exp(-h * fator * np.arange(5000)))
    assert abs(case.ajuste_h()[0] - h) / h < 1e-6


# ---------------------------------------------------------------------------
# Contador do tempo de simulacao
# ---------------------------------------------------------------------------
def test_contador_andamento(tmp_path):
    """Toda simulacao grava andamento.json (contador) e o acompanhar.py --todos o le."""
    sim, out = rodar(tmp_path / "caso", "--case", "duct", "--res", "8", "--max-iter", "300")
    a = json.load(open(os.path.join(out, "andamento.json")))
    assert a["situacao"] in ("terminou", "interrompido") and a["it"] == sim.it
    assert a["decorrido_total_s"] > 0
    assert os.path.exists(os.path.join(out, "graficos", "andamento.json"))
    log = open(os.path.join(out, "execucao.log")).read()
    assert "[contador]" in log and "tempo total de simulacao" in log
    linhas = acompanhar.tabela_contadores(acompanhar.andamentos([str(tmp_path)]), todas=True)
    assert any("duct" in l for l in linhas)
    html = open(os.path.join(out, "graficos", "convergencia.html"), encoding="utf-8").read()
    assert "contadorTxt" in html and '"andamento": {' in html


# ---------------------------------------------------------------------------
# Robustez: nunca 'list index out of range'
# ---------------------------------------------------------------------------
def test_sem_indexerror_salvar_antes_do_primeiro_relatorio(tmp_path):
    """Salvar parcial / graficos antes do primeiro relatorio (historico vazio)."""
    a = lbm3d.parse_args(["--case", "duct", "--res", "8", "--no-gui"])
    case = lbm3d.CASES["duct"](a)
    sim = lbm3d.Simulation(case, a)
    sim.step()
    pastas = lbm3d.Pastas(str(tmp_path))
    lbm3d.salvar_parcial(sim, pastas, a)
    assert figuras.graficos_convergencia(sim, str(tmp_path)) == []
    sim.hist["C_D"] = []                      # grandeza sem valores
    sim.hist["it"].append(1)
    sim.hist["res"].append(1.0)
    figuras.graficos_convergencia(sim, str(tmp_path))
    assert sim.residuos()["res"] == 1.0
    assert lbm3d._ultimo([]) != lbm3d._ultimo([])   # nan


def test_sem_indexerror_csv_com_colunas_variaveis(tmp_path):
    """Grandezas que aparecem/somem no meio da execucao nao desalinham o CSV, e o
    acompanhar.py le CSVs cortados, so com cabecalho ou com colunas a mais."""
    a = lbm3d.parse_args(["--case", "duct", "--res", "8", "--no-gui"])
    sim = lbm3d.Simulation(lbm3d.CASES["duct"](a), a)
    g = str(tmp_path)
    for n in range(3):
        sim.step()
        sim.record()
        if n == 1:
            sim.hist["nova"] = [1.0]
        lbm3d.log_convergence(g, sim, a, 1.0)
    linhas = [l for l in open(os.path.join(g, "convergencia.csv")) if not l.startswith("#")]
    assert len({len(l.split(",")) for l in linhas}) == 1
    meta, c = acompanhar.ler(g)
    assert c is not None and len(c["it"]) == 3
    with open(os.path.join(g, "convergencia.csv"), "a") as fh:
        fh.write("4,1e-3")                    # linha cortada (arquivo sendo gravado)
    meta, c = acompanhar.ler(g)
    assert len(c["it"]) == 4
    acompanhar.texto_status(meta, c, acompanhar.diagnostico(meta, c))
    so_cab = tmp_path / "vazio"
    so_cab.mkdir()
    (so_cab / "convergencia.csv").write_text("# caso=x tol=1e-7\nit,res,tempo_unix,mlups\n")
    assert acompanhar.ler(str(so_cab))[1] is None
    acompanhar.gerar_html(str(so_cab))


def test_conjugado_temperatura_uniforme_se_mantem():
    """Com solido conjugado e T = Tinf em todo o dominio, T continua uniforme mesmo com o
    escoamento partindo do repouso abrupto (u = U encostado no corpo).  Com fluxo advectivo
    T u (Tref = 0) os nos da interface desviavam ate ~80% (fonte -T div u, onde u cai a zero);
    com (T - Tinf) u (Tref = Tinf) T = Tinf e solucao exata."""
    a = lbm3d.parse_args(["--case", "calor", "--res", "8", "--Ds", "10", "--Re", "100",
                          "--metodo", "com_solido", "--no-gui"])
    case = lbm3d.CASES["calor"](a)
    sim = lbm3d.Simulation(case, a)
    th = sim.termico
    assert th.Tref == 1.0
    th.T[...] = np.where(th.ativo.astype(bool), 1.0, 0.0)
    u0 = np.where(case.solid[None], 0.0, sim.u)
    th.g = np.ascontiguousarray(lbm3d.equilibrio_termico(th.T, th.C, u0, th.Tref))
    th.g[:, ~th.ativo.astype(bool)] = 0.0
    for _ in range(400):
        sim.step()
    T = th.temperatura()
    assert np.nanmax(np.abs(T - 1.0)) < 1e-10, np.nanmax(np.abs(T - 1.0))
