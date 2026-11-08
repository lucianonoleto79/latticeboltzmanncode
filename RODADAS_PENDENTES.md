# Rodadas pendentes do cilindro (para rodar no seu computador)

Preparado em 08/11/2026. Tempos estimados a partir das rodadas com D = 6 no contêiner (~1,2 MLUPS);
no seu computador a rodada antiga do cilindro fez 0,7–1,0 MLUPS, então conte com até ~1,7 × mais.
Cada processo usa basicamente um núcleo: rodadas em terminais separados podem andar ao mesmo tempo.

## 0. Preparar o código (uma vez)

Com o bundle baixado da sessão (`latticeboltzmanncode.bundle`):

```bash
git clone latticeboltzmanncode.bundle latticeboltzmanncode
cd latticeboltzmanncode
git checkout claude/beautiful-brahmagupta-tmdwf3
pip install -r requirements.txt
```

No macOS, `caffeinate -i` na frente de cada comando impede que o computador durma durante a rodada.

## 1. Rodadas

| # | O que é | Malha | Tempo estimado (contêiner → seu computador) | Comando |
|---|---|---|---|---|
| A | Convergência de malha, D = 10 (Um = 0,05) | 251 × 43 × 43 (464 mil nós) | ~3 h → ~5 h | `caffeinate -i bash scripts/estudo_malha_cilindro.sh 10` |
| B | 3D-2Z (Re = 100, transiente), D = 10 | 251 × 43 × 43 (464 mil nós) | ~3,5 h → ~6 h | `caffeinate -i python lbm3d.py --case cylinder2z --res 10 --no-gui` |
| C | Convergência de malha, D = 16 (Um = 0,05) | 401 × 68 × 68 (1,85 milhão) | ~18 h → ~30 h | `caffeinate -i bash scripts/estudo_malha_cilindro.sh 16` |
| D | Convergência de malha, D = 20 (Um = 0,05) | 501 × 84 × 84 (3,5 milhões) | ~45 h → ~75 h | `caffeinate -i bash scripts/estudo_malha_cilindro.sh 20` |

Prioridade sugerida: A e B (respondem às perguntas principais: o C_L entra na faixa com o refino? há
desprendimento de vórtices com D = 10?), depois C; D só se C ainda não for conclusivo.

As estimativas de iterações supõem que elas crescem com D/Um (calibrado no D = 6, que convergiu ~20 %
antes do previsto) e podem errar ±30 %. B roda um número fixo de iterações (150 tempos convectivos =
33 800 it. com D = 10).

## 2. Durante a rodada

* Acompanhar: abra `casos/<pasta>/graficos/convergencia.html` no navegador (recarrega sozinho);
  o terminal mostra a tabela de resíduos a cada 100 it.
* Interromper: `Ctrl+C` abre o menu (continuar / salvar parcial / pausar / salvar e sair), ou, de
  outro terminal: `python acompanhar.py casos/<pasta> --parar` (também `--salvar`, `--pausar`,
  `--continuar`).
* Continuar depois: rode o mesmo comando de novo (o script do estudo de malha retoma sozinho do
  checkpoint e pula o que já terminou), ou, para o 3D-2Z:

```bash
caffeinate -i python lbm3d.py --restart casos/cylinder2z_D10_Re100 --no-gui
```

Pastas geradas: `casos/cylinder_D10_Re20_U0.05`, `casos/cylinder_D16_Re20_U0.05`,
`casos/cylinder_D20_Re20_U0.05`, `casos/cylinder2z_D10_Re100`.

## 3. Depois de terminar

Análise (pode rodar você mesmo; usa o D = 6 que já está em `resultados/`):

```bash
python scripts/convergencia_malha.py
python scripts/pos_cilindro.py casos/cylinder_D10_Re20_U0.05
```

Para eu processar e documentar (README, `resultados/README.md`), empacote as pastas sem checkpoint
nem `.vtk` e envie o arquivo na conversa:

```bash
tar --exclude=checkpoint --exclude='*.vtk' --exclude='parcial_*' -czf rodadas_cilindro.tar.gz casos/cylinder_D*_Re20_U0.05 casos/cylinder2z_D10_Re100
```

Se quiser que eu calcule também a bolha de recirculação (`pos_cilindro.py`) com o campo, mande à
parte o `.vtk` de `escoamento/final/` (D = 10: ~9 MB).
