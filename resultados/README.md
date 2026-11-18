# Casos simulados

Cada pasta segue a organização descrita no [README principal](../README.md#organização-dos-casos):
`dominio/`, `contorno/`, `malha/`, `escoamento/`, `graficos/`, mais `parametros.json`, `execucao.log`,
`estado.json`, `validacao.txt` e `LEIAME.txt`. Checkpoints (`.npz`) e campos `.vtk` não são versionados.
Todos os casos usam o modelo K17 e os parâmetros do benchmark correspondente.

| Pasta | Caso | Resultado principal | Reproduzir |
|---|---|---|---|
| [`duct_24x16`](duct_24x16) | Duto retangular 24 × 16, ν = 0,02 | erro L2 = 1,694·10⁻³ (White 2006) | `python lbm3d.py --case duct` |
| [`cylinder_D10_Re20`](cylinder_D10_Re20) | Cilindro Schäfer–Turek 3D-1Z, D = 10, Re = 20 | C_D = 6,230 ✔, Δp = 0,1654 Pa ✔, C_L = 0,0102 (fora da faixa) | `python lbm3d.py --case cylinder --res 10` |
| [`cylinder_D6_Re20_U0.05`](cylinder_D6_Re20_U0.05) | Cilindro 3D-1Z, D = 6, Um = 0,05 (estudo de malha) | C_D = 6,329 ✘, C_L = 0,0235 ✘, Δp = 0,1656 Pa ✔ | `python lbm3d.py --case cylinder --res 6 --u-lattice 0.05 --tol 1e-6` |
| [`cylinder_convergencia_malha`](cylinder_convergencia_malha) | Estudo de malha do cilindro (Um = 0,05) | em andamento: só D = 6 | `python scripts/convergencia_malha.py` |
| [`cylinder_D6_Re20_U0.05_cumulant`](cylinder_D6_Re20_U0.05_cumulant) | Cilindro 3D-1Z, D = 6, Um = 0,05, cumulante clássico | C_D = 6,764 ✘, C_L = 0,0337 ✘, Δp = 0,1811 Pa ✘ | `python lbm3d.py --case cylinder --res 6 --u-lattice 0.05 --tol 1e-6 --model cumulant` |
| [`cylinder_D6_Re20_U0.05_bgk`](cylinder_D6_Re20_U0.05_bgk) | Cilindro 3D-1Z, D = 6, Um = 0,05, BGK | C_D = 6,883 ✘, C_L = 0,0294 ✘, Δp = 0,1860 Pa ✘ | `python lbm3d.py --case cylinder --res 6 --u-lattice 0.05 --tol 1e-6 --model bgk` |
| [`cylinder_comparacao_modelos`](cylinder_comparacao_modelos) | Comparação K17 × cumulante × BGK (D = 6) | K17 o mais próximo nas três grandezas | `python scripts/comparar_cilindro.py ...` (ver README) |
| [`cylinder_D6_Re20`](cylinder_D6_Re20) | Cilindro 3D-1Z, D = 6, Um = 0,1 (estudo de Mach) | C_D = 6,302 ✘, C_L = 0,0214 ✘, Δp = 0,1642 Pa ✘ (0,5 % abaixo) | `python lbm3d.py --case cylinder --res 6 --tol 1e-6` |
| [`cylinder_D6_Re20_U0.025`](cylinder_D6_Re20_U0.025) | Cilindro 3D-1Z, D = 6, Um = 0,025 (estudo de Mach) | C_D = 6,343 ✘, C_L = 0,0246 ✘, Δp = 0,1666 Pa ✔ | `python lbm3d.py --case cylinder --res 6 --u-lattice 0.025 --tol 1e-6` |
| [`cylinder_efeito_mach`](cylinder_efeito_mach) | Efeito de Um (Mach) no cilindro, D = 6 | variação ≤ 0,65 % em C_D, 15 % em C_L | `python scripts/comparar_cilindro.py ...` (ver README) |
| [`cylinder_D6_Re20_U0.05_lim0.01`](cylinder_D6_Re20_U0.05_lim0.01) | Cilindro 3D-1Z, D = 6, Um = 0,05, K17 com limitador 0,01 | C_D = 6,367 ✘, C_L = 0,0241 ✘, Δp = 0,1667 Pa ✔ | `python lbm3d.py --case cylinder --res 6 --u-lattice 0.05 --tol 1e-6 --limiter 0.01` |
| [`cylinder_efeito_limitador`](cylinder_efeito_limitador) | Efeito do limitador (D = 6) | limitador piora C_D em 0,6 % | `python scripts/comparar_cilindro.py ...` (ver README) |
| [`cylinder_pos_processamento`](cylinder_pos_processamento) | L_a, separação e C_p das rodadas D = 6 | K17: L_a ≈ 0,66 D, separação ≈ 151° | `python scripts/pos_cilindro.py casos/<caso>` |
| [`cylinder2z_D6_Re100`](cylinder2z_D6_Re100) | Cilindro 3D-2Z (Re = 100, transiente), D = 6 | sem desprendimento de vórtices; C_D = 3,66, Δp = 3,40 Pa | `python lbm3d.py --case cylinder2z --res 6` |
| [`spheres_L32_chi0.3`](spheres_L32_chi0.3) | Esferas, L = 32, χ = 0,3 | C = 1,7132 (erro +0,77 %) | `python lbm3d.py --case spheres --res 32 --chi 0.3` |
| [`spheres_L32_chi0.5`](spheres_L32_chi0.5) | Esferas, L = 32, χ = 0,5 | C = 2,8583 (erro +0,47 %) | `python lbm3d.py --case spheres --res 32 --chi 0.5` |
| [`spheres_L32_chi0.7`](spheres_L32_chi0.7) | Esferas, L = 32, χ = 0,7 | C = 6,0623 (erro +1,04 %) | `python lbm3d.py --case spheres --res 32 --chi 0.7` |
| [`spheres_L48_chi0.5`](spheres_L48_chi0.5) | Esferas, L = 48, χ = 0,5 | C = 2,8516 (erro +0,23 %) | `python lbm3d.py --case spheres --res 48 --chi 0.5` |

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

## cylinder_D6_Re20_U0.05 — Cilindro 3D-1Z, D = 6, Um = 0,05

Primeira malha do estudo de convergência (K17, ω = 1,923, ν = 0,00667, malha 151 × 27 × 27). Convergiu
(resíduo 7,0·10⁻⁷ < 10⁻⁶) em 15 600 it., 25 min no contêiner a 1,16 MLUPS. C_D = 6,32927 e C_L = 0,02346
fora das faixas; Δp = 0,16565 Pa dentro. Rodado com a versão atual; parciais, `.vtk` e checkpoint não
versionados.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/pos_cilindro.png`
* `graficos/pos_cilindro.txt`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## cylinder_convergencia_malha — Estudo de malha do cilindro (Um = 0,05)

Gerado por `python scripts/convergencia_malha.py` a partir das pastas `cylinder_D*_Re20_U0.05`. Por
enquanto só tem D = 6; a extrapolação de Richardson aparece com 3 malhas.

Arquivos:

* `convergencia_malha.png`
* `convergencia_malha.txt`

## cylinder_D6_Re20_U0.05_cumulant — Cilindro 3D-1Z, D = 6, Um = 0,05, cumulante clássico [G15]

Convergiu (resíduo 5,7·10⁻⁷) em 16 200 it., 25 min no contêiner. C_D = 6,76358, C_L = 0,03374,
Δp = 0,18107 Pa, os três fora das faixas. Parciais, `.vtk` e checkpoint não versionados.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/pos_cilindro.png`
* `graficos/pos_cilindro.txt`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## cylinder_D6_Re20_U0.05_bgk — Cilindro 3D-1Z, D = 6, Um = 0,05, BGK

Convergiu (resíduo 9,3·10⁻⁷) em 15 700 it., ~25 min no contêiner; estável apesar de ω = 1,923. C_D =
6,88256, C_L = 0,02935, Δp = 0,18602 Pa, os três fora das faixas.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/pos_cilindro.png`
* `graficos/pos_cilindro.txt`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## cylinder_comparacao_modelos — Comparação de modelos de colisão (D = 6)

Gerado por `scripts/comparar_cilindro.py` a partir das três pastas acima e de `cylinder_D6_Re20_U0.05`.

Arquivos:

* `comparacao_modelos.png`
* `comparacao_modelos.txt`

## cylinder_D6_Re20 — Cilindro 3D-1Z, D = 6, Um = 0,1 (padrão)

K17, ω = 1,852, ν = 0,01333. Convergiu (resíduo 9,8·10⁻⁷) em 9 900 it., 15 min no contêiner (rodando em
paralelo com outra simulação). C_D = 6,30157, C_L = 0,02135, Δp = 0,16418 Pa.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/pos_cilindro.png`
* `graficos/pos_cilindro.txt`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## cylinder_D6_Re20_U0.025 — Cilindro 3D-1Z, D = 6, Um = 0,025

K17, ω = 1,961, ν = 0,00333. Convergiu (resíduo 8,8·10⁻⁷) em 22 300 it., 35 min no contêiner. C_D =
6,34293, C_L = 0,02463, Δp = 0,16657 Pa.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/pos_cilindro.png`
* `graficos/pos_cilindro.txt`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## cylinder_efeito_mach — Efeito da velocidade de rede (D = 6)

Gerado por `scripts/comparar_cilindro.py` a partir de `cylinder_D6_Re20`, `cylinder_D6_Re20_U0.05` e
`cylinder_D6_Re20_U0.025`.

Arquivos:

* `efeito_mach.png`
* `efeito_mach.txt`

## cylinder_D6_Re20_U0.05_lim0.01 — Cilindro 3D-1Z, D = 6, Um = 0,05, limitador 0,01

Convergiu (resíduo 7,1·10⁻⁷) em 15 600 it., 25 min no contêiner. C_D = 6,36690, C_L = 0,02410,
Δp = 0,16670 Pa.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/pos_cilindro.png`
* `graficos/pos_cilindro.txt`
* `malha/malha.txt`
* `malha/malha_corte.png`
* `malha/malha_detalhe.png`

## cylinder_efeito_limitador — Efeito do limitador de 3ª ordem (D = 6)

Gerado por `scripts/comparar_cilindro.py` a partir de `cylinder_D6_Re20_U0.05` e
`cylinder_D6_Re20_U0.05_lim0.01`.

Arquivos:

* `efeito_limitador.png`
* `efeito_limitador.txt`

## cylinder_pos_processamento — Grandezas dos campos (D = 6)

Tabela gerada a partir de `scripts/pos_cilindro.py` aplicado ao `.vtk` final das seis rodadas com D = 6;
cada uma delas também tem `graficos/pos_cilindro.png` e `graficos/pos_cilindro.txt`.

Arquivos:

* `pos_processamento.txt`

## cylinder2z_D6_Re100 — Cilindro 3D-2Z (Re = 100, não estacionário), D = 6

K17, ω = 1,969, ν = 0,00267, Um = 0,1 (rede), 20 300 it. = 150 tempos convectivos = 15,0 s, ~32 min no
contêiner. Sem desprendimento de vórtices nesta malha (C_L pico a pico 3·10⁻⁴ no último terço).
`graficos/serie_temporal.csv` tem C_D, C_L e Δp a cada 5 it.; `validacao.txt` foi refeito com a análise
corrigida, e o `execucao.log` (saída original) traz o St espúrio da versão anterior. Parciais, `.vtk` e
checkpoint não versionados.

Arquivos:

* `LEIAME.txt`
* `estado.json`
* `execucao.log`
* `parametros.json`
* `validacao.txt`
* `contorno/condicoes_contorno.txt`
* `contorno/contorno_3d.png`
* `contorno/contorno_cortes.png`
* `contorno/perfil_entrada.png`
* `dominio/dominio_3d.png`
* `dominio/dominio_cortes.png`
* `escoamento/final/componentes_pressao.png`
* `escoamento/final/linhas_corrente.png`
* `escoamento/final/velocidade_cortes.png`
* `graficos/convergencia.csv`
* `graficos/convergencia.html`
* `graficos/convergencia.png`
* `graficos/grandezas.png`
* `graficos/serie_temporal.csv`
* `graficos/serie_temporal.png`
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
