# Casos simulados

Cada pasta segue a organização descrita no [README principal](../README.md#organização-dos-casos):
`dominio/`, `contorno/`, `malha/`, `escoamento/`, `graficos/`, mais `parametros.json`, `execucao.log`,
`estado.json`, `validacao.txt` e `LEIAME.txt`. Checkpoints (`.npz`) e campos `.vtk` não são versionados.
Todos os casos usam os parâmetros do benchmark correspondente e o modelo K17 (em `calor_CN26`, K17 quando ω ≥ 1,75 e cumulante clássico abaixo disso, o padrão do caso).

| Pasta | Caso | Resultado principal | Reproduzir |
|---|---|---|---|
| [`duct_24x16`](duct_24x16) | Duto retangular 24 × 16, ν = 0,02 | erro L2 = 1,694·10⁻³ (White 2006) | `python lbm3d.py --case duct` |
| [`cylinder_D10_Re20`](cylinder_D10_Re20) | Cilindro Schäfer–Turek 3D-1Z, D = 10, Re = 20 | C_D = 6,230 ✔, Δp = 0,1654 Pa ✔, C_L = 0,0102 (fora da faixa) | `python lbm3d.py --case cylinder --res 10` |
| [`spheres_L32_chi0.3`](spheres_L32_chi0.3) | Esferas, L = 32, χ = 0,3 | C = 1,7132 (erro +0,77 %) | `python lbm3d.py --case spheres --res 32 --chi 0.3` |
| [`spheres_L32_chi0.5`](spheres_L32_chi0.5) | Esferas, L = 32, χ = 0,5 | C = 2,8583 (erro +0,47 %) | `python lbm3d.py --case spheres --res 32 --chi 0.5` |
| [`spheres_L32_chi0.7`](spheres_L32_chi0.7) | Esferas, L = 32, χ = 0,7 | C = 6,0623 (erro +1,04 %) | `python lbm3d.py --case spheres --res 32 --chi 0.7` |
| [`spheres_L48_chi0.5`](spheres_L48_chi0.5) | Esferas, L = 48, χ = 0,5 | C = 2,8516 (erro +0,23 %) | `python lbm3d.py --case spheres --res 48 --chi 0.5` |
| [`calor_CN26`](calor_CN26) | Convecção forçada sobre cilindros [CN26]: 48 casos (círculo com e sem sólido em 3 malhas, 6 seções não circulares), Re = 1–200 | Nu fino sem sólido 0,820 / 1,868 / 5,125 / 7,153 (artigo 0,798 / 1,837 / 5,109 / 7,394); comparação completa em [`comparacao_CN26.md`](calor_CN26/comparacao_CN26.md) | `python campanha_calor.py --plano completo --rodar` e `--relatorio` |

## duct_24x16 — Duto retangular 24 × 16, ν = 0,02

Resultado: erro L2 = 1,694·10⁻³ (White 2006). Reproduzir: `python lbm3d.py --case duct`.

Rodado com a versão atual; todas as subpastas preenchidas. `.vtk` e checkpoint não versionados.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/duct_validacao.png`
* `graficos/grandezas.png`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## cylinder_D10_Re20 — Cilindro Schäfer–Turek 3D-1Z, D = 10, Re = 20

Resultado: C_D = 6,230 ✔, Δp = 0,1654 Pa ✔, C_L = 0,0102 (fora da faixa). Reproduzir: `python lbm3d.py --case cylinder --res 10`.

Rodado antes da reorganização: log, validação e gráfico da janela são da rodada original; `dominio/`, `contorno/`, `malha/` e `parametros.json` regerados com `--preprocess-only`. Sem `escoamento/` (campos não guardados).

Arquivos:

* `LEIAME.txt`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `graficos/convergencia.png`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## spheres_L32_chi0.3 — Esferas, L = 32, χ = 0,3

Resultado: C = 1,7132 (erro +0,77 %). Reproduzir: `python lbm3d.py --case spheres --res 32 --chi 0.3`.

Rodado antes da reorganização: log, CSV, painel HTML e gráfico da janela são da rodada original (o painel foi regerado a partir do CSV); `dominio/`, `contorno/`, `malha/` e `parametros.json` regerados com `--preprocess-only`. Sem `escoamento/` (campos não guardados).

Arquivos:

* `LEIAME.txt`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## spheres_L32_chi0.5 — Esferas, L = 32, χ = 0,5

Resultado: C = 2,8583 (erro +0,47 %). Reproduzir: `python lbm3d.py --case spheres --res 32 --chi 0.5`.

Rodado antes da reorganização: log, CSV, painel HTML e gráfico da janela são da rodada original (o painel foi regerado a partir do CSV); `dominio/`, `contorno/`, `malha/` e `parametros.json` regerados com `--preprocess-only`. Sem `escoamento/` (campos não guardados).

Arquivos:

* `LEIAME.txt`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## spheres_L32_chi0.7 — Esferas, L = 32, χ = 0,7

Resultado: C = 6,0623 (erro +1,04 %). Reproduzir: `python lbm3d.py --case spheres --res 32 --chi 0.7`.

Rodado antes da reorganização: log, CSV, painel HTML e gráfico da janela são da rodada original (o painel foi regerado a partir do CSV); `dominio/`, `contorno/`, `malha/` e `parametros.json` regerados com `--preprocess-only`. Sem `escoamento/` (campos não guardados).

Arquivos:

* `LEIAME.txt`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## spheres_L48_chi0.5 — Esferas, L = 48, χ = 0,5

Resultado: C = 2,8516 (erro +0,23 %). Reproduzir: `python lbm3d.py --case spheres --res 48 --chi 0.5`.

Rodado antes da reorganização: log, CSV, painel HTML e gráfico da janela são da rodada original (o painel foi regerado a partir do CSV); `dominio/`, `contorno/`, `malha/` e `parametros.json` regerados com `--preprocess-only`. Sem `escoamento/` (campos não guardados).

Arquivos:

* `LEIAME.txt`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## calor_CN26 — Convecção forçada sobre cilindros, campanha de Cunha & Noleto (2026)

Resultado: tabelas no formato do artigo (Tabelas 2–6, Fig. 12, Strouhal) em `comparacao_CN26.md`. Reproduzir: `python campanha_calor.py --plano completo --rodar --processos 3` e depois `python campanha_calor.py --relatorio`; o relatório também copia os casos para cá.

Rodado com a versão atual. As 48 rodadas estão completas em `casos/`, cada uma com a organização padrão (`dominio/`, `contorno/`, `malha/`, `escoamento/final/`, `graficos/`, logs e json). Ficaram de fora os checkpoints, os `.npz` e os `.vtk`. Em `resultado_calor.json`, `Strouhal` é o St pelo período de C_L; `Strouhal_Nu` é o St pelo período do Nu, que só vale para corpos simétricos.

Arquivos na raiz:

* `comparacao_CN26.md`
* `fig12_correlacoes.png`
* `nusselt_lbm_x_artigo.png`
* `resultados.json`

Arquivos de cada caso com sólido (exemplo `casos/calor_circulo_comsolido_Re200_D20/`):

* `LEIAME.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/temperatura.png`
* `escoamento/final/velocidade_cortes.png`
* `estado.json`
* `execucao.log`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/nusselt.csv`
* `graficos/solido_ajuste_h.png`
* `graficos/solido_perfil_nusselt.png`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`
* `parametros.json`
* `resultado_calor.json`
* `validacao.txt`

Os casos sem sólido não têm `solido_*.png`. Os transientes (Re = 100 e 200) têm também:

* `graficos/nusselt_7periodos.png`
* `graficos/nusselt_fft.png`
* `graficos/nusselt_historia.png`

Pastas dos casos (48):

* [`calor_circulo_comsolido_Re100_D10`](calor_CN26/casos/calor_circulo_comsolido_Re100_D10)
* [`calor_circulo_comsolido_Re100_D14`](calor_CN26/casos/calor_circulo_comsolido_Re100_D14)
* [`calor_circulo_comsolido_Re100_D20`](calor_CN26/casos/calor_circulo_comsolido_Re100_D20)
* [`calor_circulo_comsolido_Re10_D10`](calor_CN26/casos/calor_circulo_comsolido_Re10_D10)
* [`calor_circulo_comsolido_Re10_D5`](calor_CN26/casos/calor_circulo_comsolido_Re10_D5)
* [`calor_circulo_comsolido_Re10_D7`](calor_CN26/casos/calor_circulo_comsolido_Re10_D7)
* [`calor_circulo_comsolido_Re1_D10`](calor_CN26/casos/calor_circulo_comsolido_Re1_D10)
* [`calor_circulo_comsolido_Re1_D5`](calor_CN26/casos/calor_circulo_comsolido_Re1_D5)
* [`calor_circulo_comsolido_Re1_D7`](calor_CN26/casos/calor_circulo_comsolido_Re1_D7)
* [`calor_circulo_comsolido_Re200_D10`](calor_CN26/casos/calor_circulo_comsolido_Re200_D10)
* [`calor_circulo_comsolido_Re200_D14`](calor_CN26/casos/calor_circulo_comsolido_Re200_D14)
* [`calor_circulo_comsolido_Re200_D20`](calor_CN26/casos/calor_circulo_comsolido_Re200_D20)
* [`calor_circulo_semsolido_Re100_D10`](calor_CN26/casos/calor_circulo_semsolido_Re100_D10)
* [`calor_circulo_semsolido_Re100_D14`](calor_CN26/casos/calor_circulo_semsolido_Re100_D14)
* [`calor_circulo_semsolido_Re100_D20`](calor_CN26/casos/calor_circulo_semsolido_Re100_D20)
* [`calor_circulo_semsolido_Re10_D10`](calor_CN26/casos/calor_circulo_semsolido_Re10_D10)
* [`calor_circulo_semsolido_Re10_D5`](calor_CN26/casos/calor_circulo_semsolido_Re10_D5)
* [`calor_circulo_semsolido_Re10_D7`](calor_CN26/casos/calor_circulo_semsolido_Re10_D7)
* [`calor_circulo_semsolido_Re1_D10`](calor_CN26/casos/calor_circulo_semsolido_Re1_D10)
* [`calor_circulo_semsolido_Re1_D5`](calor_CN26/casos/calor_circulo_semsolido_Re1_D5)
* [`calor_circulo_semsolido_Re1_D7`](calor_CN26/casos/calor_circulo_semsolido_Re1_D7)
* [`calor_circulo_semsolido_Re200_D10`](calor_CN26/casos/calor_circulo_semsolido_Re200_D10)
* [`calor_circulo_semsolido_Re200_D14`](calor_CN26/casos/calor_circulo_semsolido_Re200_D14)
* [`calor_circulo_semsolido_Re200_D20`](calor_CN26/casos/calor_circulo_semsolido_Re200_D20)
* [`calor_elipse_semsolido_Re100_D20`](calor_CN26/casos/calor_elipse_semsolido_Re100_D20)
* [`calor_elipse_semsolido_Re10_D10`](calor_CN26/casos/calor_elipse_semsolido_Re10_D10)
* [`calor_elipse_semsolido_Re1_D10`](calor_CN26/casos/calor_elipse_semsolido_Re1_D10)
* [`calor_elipse_semsolido_Re200_D20`](calor_CN26/casos/calor_elipse_semsolido_Re200_D20)
* [`calor_hexagono45_semsolido_Re100_D20`](calor_CN26/casos/calor_hexagono45_semsolido_Re100_D20)
* [`calor_hexagono45_semsolido_Re10_D10`](calor_CN26/casos/calor_hexagono45_semsolido_Re10_D10)
* [`calor_hexagono45_semsolido_Re1_D10`](calor_CN26/casos/calor_hexagono45_semsolido_Re1_D10)
* [`calor_hexagono45_semsolido_Re200_D20`](calor_CN26/casos/calor_hexagono45_semsolido_Re200_D20)
* [`calor_hexagono_semsolido_Re100_D20`](calor_CN26/casos/calor_hexagono_semsolido_Re100_D20)
* [`calor_hexagono_semsolido_Re10_D10`](calor_CN26/casos/calor_hexagono_semsolido_Re10_D10)
* [`calor_hexagono_semsolido_Re1_D10`](calor_CN26/casos/calor_hexagono_semsolido_Re1_D10)
* [`calor_hexagono_semsolido_Re200_D20`](calor_CN26/casos/calor_hexagono_semsolido_Re200_D20)
* [`calor_placa_semsolido_Re100_D20`](calor_CN26/casos/calor_placa_semsolido_Re100_D20)
* [`calor_placa_semsolido_Re10_D10`](calor_CN26/casos/calor_placa_semsolido_Re10_D10)
* [`calor_placa_semsolido_Re1_D10`](calor_CN26/casos/calor_placa_semsolido_Re1_D10)
* [`calor_placa_semsolido_Re200_D20`](calor_CN26/casos/calor_placa_semsolido_Re200_D20)
* [`calor_quadrado45_semsolido_Re100_D20`](calor_CN26/casos/calor_quadrado45_semsolido_Re100_D20)
* [`calor_quadrado45_semsolido_Re10_D10`](calor_CN26/casos/calor_quadrado45_semsolido_Re10_D10)
* [`calor_quadrado45_semsolido_Re1_D10`](calor_CN26/casos/calor_quadrado45_semsolido_Re1_D10)
* [`calor_quadrado45_semsolido_Re200_D20`](calor_CN26/casos/calor_quadrado45_semsolido_Re200_D20)
* [`calor_quadrado_semsolido_Re100_D20`](calor_CN26/casos/calor_quadrado_semsolido_Re100_D20)
* [`calor_quadrado_semsolido_Re10_D10`](calor_CN26/casos/calor_quadrado_semsolido_Re10_D10)
* [`calor_quadrado_semsolido_Re1_D10`](calor_CN26/casos/calor_quadrado_semsolido_Re1_D10)
* [`calor_quadrado_semsolido_Re200_D20`](calor_CN26/casos/calor_quadrado_semsolido_Re200_D20)
