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
   * Rodadas longas no contêiner só em trechos com `--max-time <min>` e `--restart`, e **sem versionar
     o checkpoint no GitHub** (com D = 10 ele já tem ~100 MB; ver "Uso do GitHub"): se o contêiner for
     recolhido entre trechos, o checkpoint se perde. Prefira o computador do usuário.
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
* Siga também as regras da seção seguinte em todo push, PR ou interação com o GitHub.

## Uso do GitHub: evitar nova suspensão da conta (vale para toda simulação e toda sessão)

Em 11/2026 a conta foi suspensa porque "some activity on your account was flagged by our abuse-detection
systems for manual review, as it may conflict with our Terms of Service" (resposta do suporte; a
atividade exata não foi informada). Regras tiradas dos termos do GitHub (Terms of Service, seções B.3
e H; Acceptable Use Policies, seções 4 e 9; termos do Actions; "About large files on GitHub"),
conferidos em 10/11/2026:

1. **Pouca atividade automática, em lotes.** Os termos proíbem "automated excessive bulk activity" e
   "undue burden on our servers through automated means".
   * Junte o trabalho de uma etapa em commits locais e faça **um push por etapa concluída**, não um por
     commit nem por arquivo.
   * Se um push ou chamada falhar por erro de rede, tente **no máximo mais uma vez**. Se a resposta for
     403, conta suspensa ou acesso negado, **pare todas as ações no GitHub** e avise o usuário. Não
     insista e não contorne.
   * Nada de laços, consultas periódicas (polling) ou check-ins agendados à API do GitHub, e nada de
     acompanhar PRs automaticamente, a menos que o usuário peça.
   * No máximo um PR por tarefa. Edite a descrição do PR uma vez por etapa, não a cada commit.
2. **Nada de interação inautêntica.** Não comente, não abra issues, não dê estrela, não siga ninguém e
   não marque pessoas (`@nome`) sem pedido explícito do usuário. Não crie contas, bots ou tokens novos e
   não use credenciais de terceiros. Os commits mantêm a atribuição ao Claude (`Co-Authored-By`), sem
   se passar por outra pessoa.
3. **GitHub Actions só para testes rápidos do código** (ex.: `pytest`). Os termos do Actions proíbem
   nos runners do GitHub "any other activity unrelated to the production, testing, deployment, or
   publication of the software project" e carga "disproportionate to the benefits". **Não rode
   simulações longas no Actions** nem fluxos que disparam a si mesmos em cadeia. Simulações longas
   rodam no computador do usuário (ver "Antes de rodar uma simulação").
4. **Tamanho dos dados** ("About large files on GitHub"): o Git avisa acima de 50 MiB e o GitHub
   bloqueia acima de 100 MiB; repositório idealmente abaixo de 1 GB.
   * Versione só o que a documentação usa: PNG, logs, CSV pequenos, `validacao.txt`, JSON.
   * Nunca versione `.npz`, `.vtk` ou checkpoints. Nenhum arquivo acima de ~10 MB sem perguntar.
   * Não versione parciais (`escoamento/parcial_*`). Dados grandes vão para fora do repositório (ou
     para um release, se o usuário pedir).
5. **Antes de todo push**, rode e confira:

   ```bash
   git diff --stat origin/<branch>..HEAD | tail -1          # quantidade de arquivos
   git diff --name-only origin/<branch>..HEAD | xargs -r ls -l 2>/dev/null | sort -k5 -n | tail -3   # maiores
   git ls-files | grep -E '\.npz$|\.vtk$|checkpoint' || echo "ok: sem dados brutos"
   git diff origin/<branch>..HEAD | grep -iE '^\+.*(api[_-]?key|secret|password|token=|BEGIN .*PRIVATE)' || echo "ok: sem segredos"
   ```

   Confira também se há caminhos absolutos de pastas temporárias ou dados pessoais em arquivos novos
   (troque por caminhos relativos).
6. Em caso de dúvida sobre se uma ação no GitHub é permitida, **pergunte ao usuário antes**.
