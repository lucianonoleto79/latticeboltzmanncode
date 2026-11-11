# Instruções para o Claude neste repositório

Solver: `lbm3d.py` (LBM D3Q27, cumulante K17). Figuras: `figuras.py`. Acompanhamento/controle:
`acompanhar.py`. Casos disponíveis: classes em `lbm3d.py` (dicionário `CASES`, opção `--case`).

## Ao simular um caso (novo ou já existente)

Toda simulação executada entra na documentação, no mesmo commit dos resultados:

1. `README.md`, seção **"Casos simulados"**: uma entrada por caso, no formato das existentes — o que
   foi simulado e com quais parâmetros, o resultado principal em uma ou duas frases (com os números),
   o comando que reproduz, o link para a pasta em `resultados/` e uma figura representativa
   (`![...](resultados/...png)`). Se o caso já tem entrada, atualize-a. Marque como ilustrativo o que
   usa parâmetros não calibrados e diga o que não foi verificado.
2. `resultados/README.md`: a tabela de casos e a seção do caso, com os arquivos gerados.
3. `README.md`, tabela **"Casos (`--case`)"**: o caso novo, com descrição e modelo padrão.

As figuras citadas no README precisam estar versionadas (o `.gitignore` libera PNG e log só dentro
de `resultados/`).

## Contador de tempo de simulação (obrigatório em toda simulação)

Toda simulação feita com este código, agora e no futuro, tem o contador de tempo presente:

* O `lbm3d.py` sempre mostra o contador (classe `Contador`): linha viva no terminal, linha
  `[contador]` a cada minuto no `execucao.log`, `andamento.json` na pasta do caso (regravado a
  cada segundo), título da janela gráfica e cabeçalho do painel `graficos/convergencia.html`.
  Não remova nem desligue o contador ao mexer no código; um caso novo herda o contador
  automaticamente (ele fica no laço principal, `_executar`).
* Ao rodar simulações (uma ou várias, em primeiro ou segundo plano), abra para o usuário um
  terminal com o contador de **todas** as simulações em curso:
  `python acompanhar.py --todos` (atualiza a cada segundo; `--raiz <pasta>` para casos fora de
  `casos/`). Ao dar notícias do andamento, informe o tempo decorrido de cada simulação em curso
  (`python acompanhar.py --todos --uma-vez`).
* Caso o contador falhe, corrija antes de continuar: a simulação não deve rodar sem ele.

## Antes de rodar uma simulação: estimar o tempo e perguntar onde rodar

Antes de qualquer simulação (caso novo ou nova rodada de um existente):

1. Estime o tempo de execução e diga como chegou nele: rode poucos passos
   (`python lbm3d.py --case <caso> ... --no-gui --max-iter 200 --outdir <scratch>`) ou use uma rodada
   anterior do mesmo caso (`execucao.log` / `graficos/convergencia.csv` em `resultados/`), meça o tempo
   por iteração (ou MLUPS) e multiplique pelo número de iterações esperado. No LBM o passo de tempo é
   `dt = 1` (unidade de rede); o número de iterações vem de `--max-iter` ou, para casos estacionários,
   da extrapolação do resíduo até `--tol` (taxa em décadas/1000 it., como no painel). Informe também o
   tamanho da malha (nx × ny × nz e número de nós).
2. Pergunte ao usuário se pode rodar, informando a estimativa. Para rodadas longas, sugira dividir
   em trechos com `--max-time <min>` e `--restart casos/<caso>` (o reinício é bit a bit idêntico a
   uma execução contínua).
3. Só comece a simulação depois da resposta.

## Casos novos

* Cada execução grava em `casos/<caso>_<parâmetros>/` (ignorada pelo git) com `dominio/`, `contorno/`,
  `malha/`, `escoamento/`, `graficos/` e `checkpoint/`; domínio, contorno e malha são gerados antes
  da simulação. A cópia versionada vai para `resultados/<caso>_<parâmetros>/`, sem `.npz` nem `.vtk`.
* Antes de uma rodada longa, gere e mostre ao usuário o domínio, o contorno e a malha com
  `python lbm3d.py --case <caso> [parâmetros] --preprocess-only`.
* Se as figuras genéricas (`figuras.py`) não bastarem, o caso define `figura_contorno(pasta, figura,
  salvar)` (figuras extras em `contorno/`, ex.: perfil de entrada do cilindro), `descricao_contorno()`
  e, se precisar, `descricao_malha()`.

## Transferência de calor (`--case calor`)

* Casos do artigo [CN26] (convecção forçada 2D sobre cilindros): `--forma`, `--Re`, `--metodo
  sem_solido|com_solido`, `--res` (D em nós), `--u-lattice`.  Os critérios de parada são os do
  artigo (Eq. 16 em janela de 5 D/U para os estacionários; duração em D/U para os transientes;
  1 − θ ≥ 0,98 no método com sólido).  O resultado vai para `resultado_calor.json`.
* O núcleo em C (`lbm_kernel.c`, compilado em `build/`) é usado automaticamente; qualquer
  mudança na colisão precisa ser feita nos dois (NumPy e C) e conferida por `tests/test_calor.py`.

## Git

* Trabalhe numa branch e faça o merge na `main` ao terminar; não faça commit direto na `main`.
* O repositório é apenas local (sem remote).
* Os dados brutos (`.npz`) e os campos `.vtk` não vão para o repositório.
