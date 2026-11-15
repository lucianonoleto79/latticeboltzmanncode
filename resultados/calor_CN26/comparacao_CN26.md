# Comparação com Cunha & Noleto (2026) [CN26]

Resultados do lbm3d.py (LBM: cumulantes D3Q27 + energia D3Q7 TRT, 2D com nz = 1) ao lado dos valores do artigo (CVFEM). Mesma disposição das tabelas do artigo; os valores do LBM estão nas colunas marcadas **LBM**. Em vez de Coarse/Medium/Fine com malhas triangulares, o LBM usa D = nós por diâmetro (malha cartesiana uniforme).

## Tabela 2 — Sensibilidade de malha, método sem sólido

| Case | Re | Mesh ID | D (nós) | Nu **LBM** | dif_i,j (%) **LBM** | Nu artigo | dif_i,j (%) artigo |
|---|---|---|---|---|---|---|---|
| 01 | 1 | Coarse | 5 | 0.819 | – | 0.731 | – |
| 01 | 1 | Medium | 7 | 0.821 | 0.25 | 0.775 | 5.68 |
| 01 | 1 | Fine | 10 | 0.820 | 0.12 | 0.798 | 2.88 |
| 02 | 10 | Coarse | 5 | 1.843 | – | 1.737 | – |
| 02 | 10 | Medium | 7 | 1.862 | 1.06 | 1.803 | 3.66 |
| 02 | 10 | Fine | 10 | 1.868 | 0.30 | 1.837 | 1.85 |
| 03 | 100 | Coarse | 10 | 5.113 | – | 4.644 | – |
| 03 | 100 | Medium | 14 | 5.104 | 0.17 | 4.987 | 6.88 |
| 03 | 100 | Fine | 20 | 5.125 | 0.41 | 5.109 | 2.39 |
| 04 | 200 | Coarse | 10 | 7.537 | – | 6.539 | – |
| 04 | 200 | Medium | 14 | 7.018 | 7.40 | 7.245 | 9.74 |
| 04 | 200 | Fine | 20 | 7.153 | 1.89 | 7.394 | 2.02 |

## Tabela 3 — Sensibilidade de malha, método com sólido

| Case | Re | Mesh ID | D (nós) | Nu **LBM** | dif_i,j (%) **LBM** | Nu artigo | dif_i,j (%) artigo |
|---|---|---|---|---|---|---|---|
| 01 | 1 | Coarse | 5 | 0.808 | – | 0.726 | – |
| 01 | 1 | Medium | 7 | 0.792 | 1.93 | 0.772 | 5.96 |
| 01 | 1 | Fine | 10 | 0.778 | 1.89 | 0.793 | 2.65 |
| 02 | 10 | Coarse | 5 | 1.871 | – | 1.815 | – |
| 02 | 10 | Medium | 7 | 1.815 | 3.13 | 1.852 | 2.00 |
| 02 | 10 | Fine | 10 | 1.756 | 3.33 | 1.863 | 0.59 |
| 03 | 100 | Coarse | 10 | 4.283 | – | 4.330 | – |
| 03 | 100 | Medium | 14 | 4.725 | 9.37 | 5.022 | 13.78 |
| 03 | 100 | Fine | 20 | 4.931 | 4.17 | 5.082 | 1.18 |
| 04 | 200 | Coarse | 10 | 5.558 | – | 5.251 | – |
| 04 | 200 | Medium | 14 | 6.667 | 16.63 | 7.266 | 27.73 |
| 04 | 200 | Fine | 20 | 7.211 | 7.55 | 7.414 | 2.00 |

## Tabela 4 — Método com sólido (malha fina)

| Case | Re | Nu_ref | Nu_num **LBM** | dif_num,ref **LBM** | Biot **LBM** | Nu_num artigo | dif artigo | Biot artigo | Sólido (LBM: (ρc_p)s/(ρc_p)f, ks/kf) |
|---|---|---|---|---|---|---|---|---|---|
| 01 | Re = 1 | 0.791 | 0.778 | 1.69% | 2.48e-03 | 0.793 | 0.25% | 2.52e-05 | 277, 157 |
| 02 | Re = 10 | 1.853 | 1.756 | 5.20% | 5.59e-03 | 1.863 | 0.54% | 5.89e-05 | 104, 157 |
| 03 | Re = 100 | 5.210 | 4.931 | 5.35% | 1.57e-02 | 5.082 | 2.46% | 1.66e-02 | 16.9, 157 |
| 04 | Re = 200 | 7.244 | 7.211 | 0.45% | 2.30e-02 | 7.414 | 2.35% | 2.31e-02 | 16.9, 157 |

## Tabela 5 — Comparação dos dois métodos (malha fina)

| Case | Re | Nu_ref | com sólido: Nu **LBM** | dif **LBM** | instantes **LBM** | sem sólido: Nu **LBM** | dif **LBM** | instantes **LBM** | com sólido: Nu artigo | dif | instantes | sem sólido: Nu artigo | dif | instantes |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 01 | Re = 1 | 0.791 | 0.778 | 1.69% | 111,660 | 0.820 | 3.67% | 49,806 | 0.793 | 0.25% | 118,682 | 0.798 | 0.88% | 8,140 |
| 02 | Re = 10 | 1.853 | 1.756 | 5.20% | 74,246 | 1.868 | 0.83% | 25,783 | 1.863 | 0.54% | 353,853 | 1.837 | 0.86% | 4,468 |
| 03 | Re = 100 | 5.210 | 4.931 | 5.35% | 94,968 | 5.125 | 1.63% | 61,200 | 5.082 | 2.46% | 217,085 | 5.109 | 1.94% | 140,000 |
| 04 | Re = 200 | 7.244 | 7.211 | 0.45% | 129,885 | 7.153 | 1.25% | 52,400 | 7.414 | 2.35% | 297,635 | 7.394 | 2.07% | 120,000 |

## Tabela 6 — Nu global para cilindros não circulares (método sem sólido)

| Re | Ellipse **LBM** | artigo | dif (%) | Vertical Plate **LBM** | artigo | dif (%) | Hexagon **LBM** | artigo | dif (%) | Hexagon (tilted 45°) **LBM** | artigo | dif (%) | Square **LBM** | artigo | dif (%) | Square (tilted 45°) **LBM** | artigo | dif (%) |
|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|---|
| 1 | 0.550 | 0.545 | 0.98 | 0.954 | 0.892 | 6.92 | 0.758 | 0.721 | 5.08 | 0.813 | 0.793 | 2.47 | 0.670 | 0.631 | 6.14 | 0.856 | 0.800 | 7.05 |
| 10 | 1.353 | 1.339 | 1.03 | 2.007 | 1.827 | 9.86 | 1.748 | 1.644 | 6.30 | 1.837 | 1.787 | 2.82 | 1.547 | 1.443 | 7.19 | 1.936 | 1.770 | 9.37 |
| 100 | 3.707 | 3.714 | 0.19 | 6.047 | 5.738 | 5.39 | 4.951 | 4.623 | 7.10 | 5.066 | 4.982 | 1.68 | 4.217 | 3.875 | 8.83 | 5.456 | 5.066 | 7.71 |
| 200 | 5.157 | 5.210 | 1.02 | 9.247 | 8.516 | 8.59 | 7.178 | 6.493 | 10.55 | 7.237 | 7.268 | 0.42 | 5.976 | 5.178 | 15.40 | 8.008 | 7.420 | 7.92 |

## Fig. 12 — Correlações Nu = C·Re^m·Pr^(1/3), 1 ≤ Re ≤ 200 (mínimos quadrados)

| Seção | C **LBM** | m **LBM** | C (Tabela 6 do artigo) | m (Tabela 6 do artigo) |
|---|---|---|---|---|
| Ellipse | 0.5950 | 0.4221 | 0.5873 | 0.4258 |
| Vertical Plate | 0.9684 | 0.4278 | 0.8981 | 0.4283 |
| Hexagon | 0.7988 | 0.4236 | 0.7628 | 0.4158 |
| Hexagon (tilted 45°) | 0.8589 | 0.4122 | 0.8335 | 0.4168 |
| Square | 0.7133 | 0.4124 | 0.6781 | 0.3998 |
| Square (tilted 45°) | 0.8983 | 0.4205 | 0.8333 | 0.4198 |

Os coeficientes "do artigo" foram ajustados aqui aos valores da Tabela 6 de [CN26] com o mesmo procedimento (mínimos quadrados de Nu = C Re^m Pr^(1/3)), porque a legenda da Fig. 12 não está no texto extraído do PDF.

## Escoamentos transientes (Re = 100 e 200, sem sólido, malha fina): Strouhal e C_D

| Seção | Re | St **LBM** | St(Nu) **LBM** | C_L rms **LBM** | C_D médio **LBM** | instantes |
|---|---|---|---|---|---|---|
| Circular | 100 | 0.1663 | 0.1657 | 0.232 | 1.346 | 61,200 |
| Circular | 200 | 0.1966 | 0.1965 | 0.480 | 1.345 | 52,400 |
| Ellipse | 100 | 0.1383 | 0.3037 | 0.007 | 1.018 | 61,200 |
| Ellipse | 200 | 0.1891 | 0.3512 | 0.109 | 0.860 | 52,400 |
| Vertical Plate | 100 | 0.1825 | 0.1828 | 0.243 | 2.139 | 61,200 |
| Vertical Plate | 200 | 0.1668 | 0.1666 | 0.345 | 2.325 | 52,400 |
| Hexagon | 100 | 0.1592 | 0.1595 | 0.252 | 1.429 | 61,200 |
| Hexagon | 200 | 0.1890 | 0.1897 | 0.456 | 1.390 | 52,400 |
| Hexagon (tilted 45°) | 100 | 0.1713 | 0.0873 | 0.228 | 1.323 | 61,200 |
| Hexagon (tilted 45°) | 200 | 0.2041 | 0.2056 | 0.466 | 1.336 | 52,400 |
| Square | 100 | 0.1502 | 0.1593 | 0.181 | 1.399 | 61,200 |
| Square | 200 | 0.1720 | 0.1727 | 0.362 | 1.343 | 52,400 |
| Square (tilted 45°) | 100 | 0.1838 | 0.1840 | 0.358 | 1.505 | 61,200 |
| Square (tilted 45°) | 200 | 0.2209 | 0.2210 | 0.660 | 1.615 | 52,400 |

St pelo período de C_L (um ciclo de desprendimento); St(Nu) = D/(2 T_Nu U), que só vale para corpos simétricos em relação ao escoamento (Nu oscila duas vezes por ciclo). No hexágono a 45° Nu oscila uma vez por ciclo. Com C_L rms ≈ 0 (elipse, Re = 100) o desprendimento é praticamente nulo e o St não é significativo.


