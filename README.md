# lbm3d — LBM 3D para escoamentos laminares (D3Q27, cumulante K17)

Solver de lattice-Boltzmann tridimensional em Python/NumPy, com operador de colisão de
cumulantes parametrizado (K17, difusão de 4ª ordem), paredes curvas, painel HTML de
convergência, checkpoint/reinício e casos de validação da literatura. Roda no seu computador
ou **na nuvem, pelo GitHub Actions** (seção abaixo).

## Estrutura do repositório

```
lbm3d.py                    solver (método, casos de teste, checkpoint, VTK, linha de comando)
acompanhar.py               acompanhamento da convergência (painel HTML, janela, terminal)
painel_template.html        modelo do painel HTML (autocontido, sem internet)
scripts/resumo_execucao.py  resumo de uma execução na nuvem (página do GitHub Actions)
tests/test_lbm3d.py         testes rápidos (pytest, ~20 s)
resultados/                 resultados de validação já obtidos (logs, CSV, gráficos, painéis)
.github/workflows/
    testes.yml              roda os testes a cada push / pull request
    simulacao.yml           roda uma simulação na nuvem, com continuação automática
requirements.txt            numpy, matplotlib
requirements-dev.txt        + pytest
```

As pastas de saída (`saida_*`, com checkpoints `.npz` e campos `.vtk`, que chegam a centenas de MB)
não vão para o repositório (`.gitignore`).

## Rodar no seu computador

Requisitos: Python ≥ 3.10, `numpy`, `matplotlib`.

```bash
pip install -r requirements.txt
python lbm3d.py --case duct                 # validação rápida (solução analítica)
python lbm3d.py --case cylinder --res 10    # benchmark Schäfer–Turek 3D-1Z
python lbm3d.py --restart saida_cylinder/checkpoint.npz   # continua de onde parou
```

Testes:

```bash
pip install -r requirements-dev.txt
pytest
```

## Rodar na nuvem (GitHub Actions)

1. Envie o repositório para o GitHub.
2. Na aba **Actions**, escolha **Simulacao LBM → Run workflow** e preencha o formulário:

   | Campo | Significado |
   |---|---|
   | `caso`, `modelo`, `res`, `chi` | como na linha de comando (`--case`, `--model`, `--res`, `--chi`) |
   | `extras` | outros argumentos, ex.: `--nu 0.01 --Re 20` |
   | `tol`, `max_iter` | critério de parada e limite total de iterações |
   | `tempo_max` | minutos de simulação por execução (padrão 320; o GitHub encerra jobs com 6 h) |
   | `reiniciar_de` | ID de uma execução anterior, para continuar do checkpoint dela |
   | `continuar` | se o tempo acabar antes de convergir, dispara sozinho a continuação |

   Também pela linha de comando, com o [GitHub CLI](https://cli.github.com):
   ```bash
   gh workflow run simulacao.yml -f caso=cylinder -f res=10 -f tol=1e-6
   ```
3. A página da execução mostra um resumo (iteração, situação, resultado da validação). Os arquivos
   ficam no artefato **`resultado`** por 30 dias: `convergencia.html` (abra no navegador),
   `convergencia.csv`/`.png`, `validacao.txt`, o último `.vtk` (ParaView), `checkpoint.npz` e o log.

**Simulações longas.** Ao atingir `tempo_max`, o solver para de forma limpa (`--max-time`), salva o
checkpoint e o artefato. Com `continuar` marcado, o fluxo dispara uma nova execução com
`reiniciar_de` apontando para a anterior; a cadeia segue até convergir ou atingir `max_iter`. O
reinício é bit a bit idêntico a uma execução contínua, e o `convergencia.csv` acumula todas as
partes. Para interromper uma cadeia, cancele a execução em andamento na aba Actions (o
solver trata o sinal de cancelamento e normalmente ainda salva o checkpoint) ou desligue o fluxo.

Os executores gratuitos (`ubuntu-latest`, 4 vCPUs) são ilimitados em repositórios públicos; em
repositórios privados contam na cota mensal de minutos da conta. A velocidade é parecida com a de
um notebook (~0,7–1 MLUPS): o cilindro com D = 10 leva ~3 h e a esfera com L = 48 ~4 h, cada um numa
execução só; casos maiores usam a continuação automática.

## Resultados de validação já obtidos

Em [`resultados/`](resultados), com o modelo K17:

| Caso | Resultado | Referência |
|---|---|---|
| Duto 24×16 (ν = 0,02) | erro L2 = 1,7·10⁻³ (BGK: 1,4·10⁻², cumulante clássico: 7,8·10⁻³) | White (2006) |
| Cilindro 3D-1Z, D = 10 | C_D = 6,230 ✔, Δp = 0,1654 Pa ✔, C_L = 0,0102 (1,7 % acima da faixa) | Schäfer & Turek (1996) |
| Esferas, L = 32, χ = 0,3 / 0,5 / 0,7 | C = 1,713 / 2,858 / 6,062 (erro +0,8 / +0,5 / +1,0 %) | Sangani & Acrivos (1982) |
| Esferas, L = 48, χ = 0,5 | C = 2,852 (erro +0,23 %; ordem observada ≈ 1,8) | idem |

## Controles durante a execução

| Ação | Janela gráfica | Sem janela (`--no-gui`) |
|---|---|---|
| Curva de convergência a cada 100 it. | gráfico à esquerda (+ fatia de \|u\|) | impressa no terminal; `convergencia.png` no fim |
| Pausar / continuar | botão **Pausar** ou tecla `P` | — |
| Salvar resultado parcial agora | botão **Salvar** ou tecla `S` | — |
| Salvar a cada 100 it. | caixa **Auto-salvar** (padrão: ligada) | `--save-every 100` (padrão) |
| Parar e reiniciar depois | botão **Salvar e sair** / `Q` / fechar janela | `Ctrl+C` |
| Reiniciar | `--restart <pasta>/checkpoint.npz` | idem |

Cada salvamento grava `<caso>_itNNNNNNN.vtk` (densidade, velocidade e máscara de sólido, para abrir
no ParaView) e sobrescreve `checkpoint.npz` (populações, iteração e histórico). O reinício reproduz
bit a bit uma execução contínua. No reinício você pode alterar `--max-iter`, `--tol`, `--save-every`,
`--no-gui`, `--outdir`, `--max-time`; os parâmetros físicos vêm do checkpoint.

`--max-time M` encerra após M minutos salvando o checkpoint (usado na nuvem); `SIGTERM` tem o mesmo
efeito que `Ctrl+C`. Ao final, `<saída>/estado.json` informa a iteração e se a simulação convergiu.

### Painel HTML de convergência
A cada 100 iterações o código acrescenta uma linha a `<saída>/convergencia.csv` e regrava
**`<saída>/convergencia.html`**, uma página autocontida (sem internet, sem bibliotecas externas).
Basta abri-la no navegador: ela se recarrega sozinha a cada 5 s. Para atualização ao vivo sem
recarregar (e para acompanhar de outro computador via túnel SSH), sirva a pasta:

```bash
python acompanhar.py saida_spheres_L32_chi0.5 --servir        # http://127.0.0.1:8765/convergencia.html
python acompanhar.py saida_spheres_L32_chi0.5                 # (re)gera o HTML e abre no navegador
```

O painel mostra:
* estado da simulação (rodando / convergida / parada ou pausada / divergiu) e há quanto tempo veio a última leitura;
* iteração, resíduo, taxa de convergência (décadas por 1000 it.), velocidade (it/s e MLUPS) e
  **estimativa do tempo restante** até a tolerância (extrapolação log-linear das últimas 20 leituras);
* curva do resíduo em escala log com a tolerância e a extrapolação;
* evolução de cada grandeza do caso (C_D, C_L, Δp, C, erro L2 …) com a **faixa de referência da
  literatura** sombreada, a variação relativa nas últimas 500 it., e opção "Tudo"/"Último terço";
* valores ao passar o mouse e uma tabela com as últimas 30 leituras; tema claro/escuro automático.

Outros modos: `--janela` (matplotlib), `--texto` (terminal, ex.: via ssh), `--png` (imagem estática).

Resíduo de convergência: `r = ||u^n − u^(n−100)||₂ / ||u^n||₂` sobre os nós de fluido; o caso para
quando `r < --tol` (padrão 1e-7).

## Método e equações

### 1. Equação de lattice-Boltzmann, rede D3Q27
`f_i(x + c_i, t+1) = f_i*(x, t)`, com `c_i ∈ {−1,0,1}³`, pesos `w = w_a w_b w_c`, sendo `(1/6, 2/3, 1/6)`,
e `c_s² = 1/3`, `ν = c_s²(1/ω₁ − 1/2)`, `p = c_s² ρ`.
Macroscópicas: `ρ = Σ f_i`, `ρu = Σ f_i c_i + F/2`.

### 2. Colisão de cumulantes [G15]
Transformação "chimera" para momentos centrais, direção por direção (z, y, x):

```
κ₀ = f₋ + f₀ + f₊
κ₁ = (f₊ − f₋) − u κ₀
κ₂ = (f₊ + f₋) − 2u (f₊ − f₋) + u² κ₀
```

Cumulantes (k = momentos centrais/ρ): `c₂ = k₂`, `c₃ = k₃`,
`c₂₁₁ = k₂₁₁ − (k₂₀₀k₀₁₁ + 2k₁₁₀k₁₀₁)`, `c₂₂₀ = k₂₂₀ − (k₂₀₀k₀₂₀ + 2k₁₁₀²)`,
`c₂₂₁ = k₂₂₁ − (k₂₀₀k₀₂₁ + k₀₂₀k₂₀₁ + 4k₁₁₀k₁₁₁ + 2k₁₀₁k₁₂₀ + 2k₀₁₁k₂₁₀)`,
`c₂₂₂ = k₂₂₂ − Σk₄k₂ − Σk₃k₃ + 2Σk₂k₂k₂` (as somas de partições estão escritas no código).

Relaxação de 2ª ordem com correção galileana [G15, Eqs. 58–63]:

```
∂x ux = −ω₁/2 (2c₂₀₀ − c₀₂₀ − c₀₀₂) − ω₂/2 (c₂₀₀ + c₀₂₀ + c₀₀₂ − 1)
(c₂₀₀ − c₀₂₀)* = (1−ω₁)(c₂₀₀ − c₀₂₀) − 3(1 − ω₁/2)(ux² ∂x ux − uy² ∂y uy)
(c₂₀₀ + c₀₂₀ + c₀₀₂)* = traço + ω₂(1 − traço) − 3(1 − ω₂/2)(ux² ∂x ux + uy² ∂y uy + uz² ∂z uz)
c₁₁₀* = (1 − ω₁) c₁₁₀
```

Força de corpo: o momento central de 1ª ordem troca de sinal (`−F/2ρ → +F/2ρ`).

### 3. Parametrização K17 (difusão de 4ª ordem) [G17], com ω₂ = 1

```
ω₃ = 8(ω₁−2)(1+2ω₁) / (7ω₁² − 14ω₁ − 8)          (c₁₂₀ + c₁₀₂, …)
ω₄ = 8(ω₁−2)(4ω₁−7) / (9ω₁² − 50ω₁ + 56)          (c₁₂₀ − c₁₀₂, …)
ω₅ = 24(ω₁−2)(3ω₁² − 7ω₁ − 2) / (29ω₁³ − 130ω₁² + 152ω₁ + 48)   (c₁₁₁)
A  = (4 + 2ω₁ − 3ω₁²) / (2 − 7ω₁ + 5ω₁²)
B  = (4 + 28ω₁ − 14ω₁²) / (6 − 21ω₁ + 15ω₁²)
c₀₂₂* = −(1/ω₁ − 1/2)·(2/3)·A·(∂y uy + ∂z uz)      (idem c₂₀₂, c₂₂₀)
c₁₁₂* = −(1/ω₁ − 1/2)·(1/3)·B·(∂x uy + ∂y ux)      (idem c₁₂₁, c₂₁₁)
c₅ = c₆ = 0
```

Essas expressões foram conferidas com o código-fonte aberto do VirtualFluids (TU Braunschweig),
de onde vêm o K17 e o artigo [G17].

**Faixa de uso:** com ω₂ = 1, A tem um polo em ω₁ = 1 e ω₄ tem um polo em ω₁ ≈ 1,556. Nos testes, o
K17 só foi estável e mais preciso para **ω₁ ≥ 1,75 (ν ≤ 0,024)**, o regime de baixa viscosidade para o
qual foi criado. Abaixo disso o código usa automaticamente o cumulante clássico [G15] e avisa.

### 4. Condições de contorno
* Parede (meio caminho) e parede móvel [L94]: `f_d̄(x_f, t+1) = f_d*(x_f) − 6 w_d ρ_w (c_d·u_w)`
* Parede curva, interpolação linear [B01], com `q` a fração do link até a parede:
  `q<½: f_d̄ = 2q f_d*(x_f) + (1−2q) f_d*(x_f − c_d)`; `q≥½: f_d̄ = f_d*/(2q) + (2q−1)/(2q) f_d̄*(x_f)`
* Saída com `ρ = 1` por extrapolação de não-equilíbrio [G02]:
  `f(x_b) = f^eq(1, u(x_{b−1})) + f(x_{b−1}) − f^eq(ρ(x_{b−1}), u(x_{b−1}))`
* Força no corpo por troca de quantidade de movimento [M02]: `F = Σ c_d [f_d*(x_f) + f_d̄(x_f, t+1)]`

## Casos de validação

| Caso | Referência | Grandeza |
|---|---|---|
| `shearwave` | decaimento exato `exp(−ν k² t)` [G17] | erro na viscosidade |
| `duct` | White (2006), *Viscous Fluid Flow*, Eq. 3-48 | erro L2 do perfil |
| `cylinder` | Schäfer & Turek (1996), caso 3D-1Z, Re = 20 | C_D ∈ [6,05; 6,25], C_L ∈ [0,008; 0,010], Δp ∈ [0,165; 0,175] Pa |
| `spheres` | Sangani & Acrivos (1982), arranjo cúbico simples, regime de Stokes | arraste C = f_t/(3πμ d u_s) |
| `cavity` | demonstração (comparar com Albensoeder & Kuhlmann 2005) | perfis centrais |

Solução do duto (−a ≤ y ≤ a, −b ≤ z ≤ b, G = −dp/dx):

```
u(y,z) = 16a²G/(μπ³) Σ_{i=1,3,5,…} (−1)^((i−1)/2) [1 − cosh(iπz/2a)/cosh(iπb/2a)] cos(iπy/2a) / i³
```

Schäfer–Turek 3D-1Z: canal 2,5 × 0,41 × 0,41 m, cilindro com D = 0,1 m em (0,5; 0,2), ν = 10⁻³ m²/s,
entrada `U = 16·Um·y·z·(H−y)(H−z)/H⁴`, com Um = 0,45 m/s, Ū = 4Um/9 = 0,2 m/s,
`C_D = 2F_x/(ρŪ²DH)`, `C_L = 2F_y/(ρŪ²DH)`, `Δp = p(0,45; 0,2; 0,205) − p(0,55; 0,2; 0,205)`.

Arranjo de esferas: célula periódica L³ com uma esfera de diâmetro d = χL (fração sólida
φ = πχ³/6), forçada por aceleração g; u_s = velocidade superficial, f_t = força de arraste por
troca de quantidade de movimento + empuxo ρgV_p. Série de Sangani & Acrivos:
`1/C = 1 − 1,7601x + x³ − 1,5593x⁶ + 3,9799x⁸ − 3,0734x¹⁰ + …`, com x = φ^(1/3) (boa até χ ≈ 0,7);
valores de referência para χ = 0,2…0,9 retirados da Tab. 1 de Bogner et al. (2015).

## Referências
* [G15] M. Geier, M. Schönherr, A. Pasquali, M. Krafczyk, *The cumulant lattice Boltzmann equation in three dimensions: Theory and validation*, Comput. Math. Appl. 70 (2015) 507–547.
* [G17] M. Geier, A. Pasquali, M. Schönherr, *Parametrization of the cumulant lattice Boltzmann method for fourth order accurate diffusion. Part I*, J. Comput. Phys. 348 (2017) 862–888.
* [L94] A.J.C. Ladd, J. Fluid Mech. 271 (1994) 285–309.
* [B01] M. Bouzidi, M. Firdaouss, P. Lallemand, Phys. Fluids 13 (2001) 3452–3459.
* [G02] Z. Guo, C. Zheng, B. Shi, Chinese Physics 11 (2002) 366–374.
* [M02] R. Mei, D. Yu, W. Shyy, L.-S. Luo, Phys. Rev. E 65 (2002) 041203.
* M. Schäfer, S. Turek, *Benchmark computations of laminar flow around a cylinder*, Notes Numer. Fluid Mech. 52 (1996) 547–566.
* A.S. Sangani, A. Acrivos, *Slow flow through a periodic array of spheres*, Int. J. Multiphase Flow 8 (1982) 343–360.
* S. Bogner, S. Mohanty, U. Rüde, *Drag correlation for dilute and moderately dense fluid-particle systems using the lattice Boltzmann method*, Int. J. Multiphase Flow 68 (2015) 71–79 ([arXiv:1401.2025](https://arxiv.org/abs/1401.2025)).
* F.M. White, *Viscous Fluid Flow*, 3ª ed., McGraw-Hill, 2006.
