# Instruções para o Claude neste repositório

Solver: `lbm3d.py` (LBM D3Q27, cumulante K17). Figuras: `figuras.py`. Acompanhamento/controle:
`acompanhar.py`. Casos disponíveis: classes em `lbm3d.py` (dicionário `CASES`, opção `--case`).

## Ao simular um caso (novo ou já existente)

Toda simulação executada entra na documentação, no mesmo commit (ou PR) dos resultados:

1. `README.md`, seção **"Casos simulados"**: uma entrada por caso, no formato das existentes — o que
   foi simulado e com quais parâmetros, o resultado principal em uma ou duas frases (com os números),
   o comando que reproduz, o link para a pasta em `resultados/` e uma figura representativa
   (`![...](resultados/...png)`). Se o caso já tem entrada, atualize-a. Marque como ilustrativo o que
   usa parâmetros não calibrados e diga o que não foi verificado.
2. `resultados/README.md`: a tabela de casos e a seção do caso, com os arquivos gerados.
3. `README.md`, tabela **"Casos (`--case`)"**: o caso novo, com descrição e modelo padrão.

As figuras citadas no README precisam estar versionadas (o `.gitignore` libera PNG e log só dentro
de `resultados/`).

## Antes de rodar uma simulação: estimar o tempo e perguntar onde rodar

Antes de qualquer simulação (caso novo ou nova rodada de um existente):

1. Estime o tempo de execução e diga como chegou nele: rode poucos passos
   (`python lbm3d.py --case <caso> ... --no-gui --max-iter 200 --outdir <scratch>`) ou use uma rodada
   anterior do mesmo caso (`execucao.log` / `graficos/convergencia.csv` em `resultados/`), meça o tempo
   por iteração (ou MLUPS) e multiplique pelo número de iterações esperado. No LBM o passo de tempo é
   `dt = 1` (unidade de rede); o número de iterações vem de `--max-iter` ou, para casos estacionários,
   da extrapolação do resíduo até `--tol` (taxa em décadas/1000 it., como no painel). Informe também o
   tamanho da malha (nx × ny × nz e número de nós).
2. Pergunte ao usuário onde rodar, recomendando conforme o tempo:
   * até ~30 min: no contêiner da sessão (os resultados são processados e commitados aqui);
   * mais de ~30 min: no computador do usuário, clonando o repositório, porque o contêiner é recolhido
     após inatividade, as tarefas em segundo plano têm limite de ~2 h e o checkpoint
     (`casos/<caso>/checkpoint/checkpoint.npz`) se perde com ele. Dê os comandos (`git clone`, o comando
     do caso e, se interromper, `python lbm3d.py --restart casos/<caso>`).
   * Alternativa para rodadas longas no contêiner: dividir em trechos com `--max-time <min>` e
     `--restart`, versionando o checkpoint numa branch entre eles.
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

## Git

* Trabalhe numa branch e abra PR para a `main`; não faça commit direto na `main`.
* Os dados brutos (`.npz`) e os campos `.vtk` não vão para o repositório.
