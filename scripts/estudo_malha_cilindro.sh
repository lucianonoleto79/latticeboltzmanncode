#!/usr/bin/env bash
# Estudo de convergencia de malha do cilindro (Schaefer & Turek 3D-1Z, Re = 20) com o K17.
# Todas as malhas usam a mesma velocidade de rede Um = 0.05 (mesmo numero de Mach; K17 valido
# ate D = 21) e o mesmo criterio de parada (residuo < 1e-6).
#
# Uso (na raiz do repositorio):
#   bash scripts/estudo_malha_cilindro.sh            # roda D = 10, 16 e 20, um depois do outro
#   bash scripts/estudo_malha_cilindro.sh 16         # so D = 16 (rode varios em terminais
#                                                    #   separados para usar mais nucleos)
# No macOS, para o computador nao dormir:  caffeinate -i bash scripts/estudo_malha_cilindro.sh 16
#
# Interromper: Ctrl+C abre o menu (salvar parcial / pausar / salvar e sair), ou, de outro terminal,
#   python acompanhar.py casos/cylinder_D16_Re20_U0.05 --parar
# Rodar o mesmo comando de novo continua do checkpoint; casos ja terminados sao pulados.
# Acompanhar: abra casos/cylinder_D16_Re20_U0.05/graficos/convergencia.html no navegador.
# Analise ao final:  python scripts/convergencia_malha.py
set -euo pipefail
cd "$(dirname "$0")/.."
[ $# -eq 0 ] && set -- 10 16 20
for D in "$@"; do
  pasta="casos/cylinder_D${D}_Re20_U0.05"
  if [ -f "$pasta/estado.json" ] && python -c "import json,sys; sys.exit(0 if json.load(open('$pasta/estado.json'))['terminou'] else 1)"; then
    echo "== D = $D ja terminou ($pasta): pulando"
  elif [ -f "$pasta/checkpoint/checkpoint.npz" ]; then
    echo "== D = $D: continuando de $pasta"
    python lbm3d.py --restart "$pasta" --no-gui
  else
    echo "== D = $D: iniciando em $pasta"
    python lbm3d.py --case cylinder --res "$D" --u-lattice 0.05 --tol 1e-6 --no-gui --save-every 5000
  fi
done
