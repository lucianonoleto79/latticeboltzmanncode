/*
 * lbm_kernel.c - nucleo em C do lbm3d.py (colisao + propagacao fundidas).
 *
 * Mesma matematica das funcoes collide/_cumulant_relax/stream do lbm3d.py (cumulante
 * classico [G15], K17 [G17], BGK em momentos centrais) e do passo termico TRT D3Q7.
 * E compilado automaticamente pelo lbm3d.py na primeira execucao (cc -O3) e usado por
 * ctypes; se nao houver compilador, o lbm3d.py usa a versao em NumPy (mesmo resultado,
 * ~10x mais lenta).
 *
 * Armazenamento: f[d*N + n], d = 9a + 3b + c (velocidade (a-1, b-1, c-1)),
 *                n = (x*ny + y)*nz + z, N = nx*ny*nz.  Propagacao periodica ("push").
 */
#include <math.h>
#include <stdint.h>

#define K(a, b, c) k[(a) * 9 + (b) * 3 + (c)]
#define BLK 512
#define BLKPAD (BLK + 24)  /* folga: evita que as 27 linhas caiam no mesmo conjunto da cache */

static const double W1[3] = {1.0 / 6.0, 2.0 / 3.0, 1.0 / 6.0};
static const double E1[3] = {1.0, 0.0, 1.0 / 3.0};

/* transformacao chimera numa linha (fm, f0, fp) -> (k0, k1, k2) */
static inline void fwd3(double *fm, double *f0, double *fp, double u)
{
    double k0 = *fm + *f0 + *fp;
    double s = *fp - *fm;
    double k2 = (*fp + *fm) - 2.0 * u * s + u * u * k0;
    double k1 = s - u * k0;
    *fm = k0; *f0 = k1; *fp = k2;
}

static inline void bwd3(double *k0, double *k1, double *k2, double u)
{
    double u2 = u * u;
    double a = *k0, b = *k1, c = *k2;
    double fm = 0.5 * ((u2 - u) * a + (2.0 * u - 1.0) * b + c);
    double fp = 0.5 * ((u2 + u) * a + (2.0 * u + 1.0) * b + c);
    double f0 = (1.0 - u2) * a - 2.0 * u * b - c;
    *k0 = fm; *k1 = f0; *k2 = fp;
}

/* eixo 0 = a, 1 = b, 2 = c */
static inline void fwd_axis(double *k, double u, int axis)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            if (axis == 0) fwd3(&K(0, i, j), &K(1, i, j), &K(2, i, j), u);
            else if (axis == 1) fwd3(&K(i, 0, j), &K(i, 1, j), &K(i, 2, j), u);
            else fwd3(&K(i, j, 0), &K(i, j, 1), &K(i, j, 2), u);
        }
}

static inline void bwd_axis(double *k, double u, int axis)
{
    for (int i = 0; i < 3; i++)
        for (int j = 0; j < 3; j++) {
            if (axis == 0) bwd3(&K(0, i, j), &K(1, i, j), &K(2, i, j), u);
            else if (axis == 1) bwd3(&K(i, 0, j), &K(i, 1, j), &K(i, 2, j), u);
            else bwd3(&K(i, j, 0), &K(i, j, 1), &K(i, j, 2), u);
        }
}

static inline double relax(double x, double w, int uselim, double lim)
{
    if (uselim && w != 1.0) {
        double ax = fabs(x);
        w = w + (1.0 - w) * ax / (lim + ax);
    }
    return (1.0 - w) * x;
}

/* prm: [model(0=bgk,1=cumulant,2=k17), w1, w2, w3, w4, w5, A, B, uselim, lim] */
static inline __attribute__((always_inline)) void cumulant_relax(double *k, double ux, double uy, double uz, const double *prm)
{
    int model = (int)prm[0];
    double w1 = prm[1], w2 = prm[2], w3 = prm[3], w4 = prm[4], w5 = prm[5];
    int uselim = (int)prm[8];
    double lim = prm[9];

    double k200 = K(2, 0, 0), k020 = K(0, 2, 0), k002 = K(0, 0, 2);
    double k110 = K(1, 1, 0), k101 = K(1, 0, 1), k011 = K(0, 1, 1);

    double T = k200 + k020 + k002;
    double dxux = -0.5 * w1 * (2.0 * k200 - k020 - k002) - 0.5 * w2 * (T - 1.0);
    double dyuy = dxux + 1.5 * w1 * (k200 - k020);
    double dzuz = dxux + 1.5 * w1 * (k200 - k002);
    double dxy = -3.0 * w1 * k110;
    double dxz = -3.0 * w1 * k101;
    double dyz = -3.0 * w1 * k011;

    double ux2 = ux * ux, uy2 = uy * uy, uz2 = uz * uz;
    double mxxMyy = (1.0 - w1) * (k200 - k020) - 3.0 * (1.0 - 0.5 * w1) * (ux2 * dxux - uy2 * dyuy);
    double mxxMzz = (1.0 - w1) * (k200 - k002) - 3.0 * (1.0 - 0.5 * w1) * (ux2 * dxux - uz2 * dzuz);
    T = T + w2 * (1.0 - T) - 3.0 * (1.0 - 0.5 * w2) * (ux2 * dxux + uy2 * dyuy + uz2 * dzuz);
    double p200 = (T + mxxMyy + mxxMzz) / 3.0;
    double p020 = (T - 2.0 * mxxMyy + mxxMzz) / 3.0;
    double p002 = (T + mxxMyy - 2.0 * mxxMzz) / 3.0;
    double p110 = (1.0 - w1) * k110;
    double p101 = (1.0 - w1) * k101;
    double p011 = (1.0 - w1) * k011;

#define PAIR(A_, B_, PA, PB) do { double P_ = relax((A_) + (B_), w3, uselim, lim); \
        double M_ = relax((A_) - (B_), w4, uselim, lim); PA = 0.5 * (P_ + M_); PB = 0.5 * (P_ - M_); } while (0)
    double p210, p012, p201, p021, p120, p102;
    PAIR(K(2, 1, 0), K(0, 1, 2), p210, p012);
    PAIR(K(2, 0, 1), K(0, 2, 1), p201, p021);
    PAIR(K(1, 2, 0), K(1, 0, 2), p120, p102);
    double p111 = relax(K(1, 1, 1), w5, uselim, lim);

    double c022 = 0, c202 = 0, c220 = 0, c112 = 0, c121 = 0, c211 = 0;
    if (model == 2) {
        double fac = -(1.0 / w1 - 0.5), A = prm[6], B = prm[7];
        c022 = fac * (2.0 / 3.0) * A * (dyuy + dzuz);
        c202 = fac * (2.0 / 3.0) * A * (dxux + dzuz);
        c220 = fac * (2.0 / 3.0) * A * (dxux + dyuy);
        c112 = fac * (1.0 / 3.0) * B * dxy;
        c121 = fac * (1.0 / 3.0) * B * dxz;
        c211 = fac * (1.0 / 3.0) * B * dyz;
    }

    double p211 = c211 + p200 * p011 + 2.0 * p110 * p101;
    double p121 = c121 + p020 * p101 + 2.0 * p110 * p011;
    double p112 = c112 + p002 * p110 + 2.0 * p101 * p011;
    double p220 = c220 + p200 * p020 + 2.0 * p110 * p110;
    double p202 = c202 + p200 * p002 + 2.0 * p101 * p101;
    double p022 = c022 + p020 * p002 + 2.0 * p011 * p011;

    double p221 = (p200 * p021 + p020 * p201 + 4.0 * p110 * p111 + 2.0 * p101 * p120 + 2.0 * p011 * p210);
    double p212 = (p200 * p012 + p002 * p210 + 4.0 * p101 * p111 + 2.0 * p110 * p102 + 2.0 * p011 * p201);
    double p122 = (p002 * p120 + p020 * p102 + 4.0 * p011 * p111 + 2.0 * p101 * p021 + 2.0 * p110 * p012);

    double S4 = (p200 * p022 + p020 * p202 + p002 * p220 + 4.0 * (p110 * p112 + p101 * p121 + p011 * p211));
    double S33 = 2.0 * (p210 * p012 + p201 * p021 + p120 * p102) + 4.0 * p111 * p111;
    double S222 = (p200 * p020 * p002 + 2.0 * (p200 * p011 * p011 + p020 * p101 * p101 + p002 * p110 * p110)
                   + 8.0 * p110 * p101 * p011);
    double p222 = S4 + S33 - 2.0 * S222;

    K(1, 0, 0) = -K(1, 0, 0);
    K(0, 1, 0) = -K(0, 1, 0);
    K(0, 0, 1) = -K(0, 0, 1);
    K(0, 0, 0) = 1.0;
    K(2, 0, 0) = p200; K(0, 2, 0) = p020; K(0, 0, 2) = p002;
    K(1, 1, 0) = p110; K(1, 0, 1) = p101; K(0, 1, 1) = p011;
    K(2, 1, 0) = p210; K(0, 1, 2) = p012; K(2, 0, 1) = p201; K(0, 2, 1) = p021;
    K(1, 2, 0) = p120; K(1, 0, 2) = p102; K(1, 1, 1) = p111;
    K(2, 2, 0) = p220; K(2, 0, 2) = p202; K(0, 2, 2) = p022;
    K(2, 1, 1) = p211; K(1, 2, 1) = p121; K(1, 1, 2) = p112;
    K(2, 2, 1) = p221; K(2, 1, 2) = p212; K(1, 2, 2) = p122;
    K(2, 2, 2) = p222;
}

/* Colisao de um bloco de nb nos (buf[d][i]); laco sem desvios por no, para o
 * compilador vetorizar (o modelo e constante dentro do laco). */
static inline __attribute__((always_inline))
void bloco(double (*buf)[BLKPAD], const uint8_t *sb, int nb, double Fx, double Fy, double Fz,
           const double *prm, const double *Wd, const double *Keq,
           double *ro, double *xo, double *yo, double *zo, const int model)
{
    const double w1 = prm[1];
#pragma clang loop vectorize(enable) interleave(enable)
    for (int i = 0; i < nb; i++) {
        double k[27];
        double rho = 0, jx = 0, jy = 0, jz = 0;
        for (int d = 0; d < 27; d++) {
            double v = buf[d][i];
            k[d] = v;
            rho += v;
        }
        for (int b = 0; b < 3; b++)
            for (int c = 0; c < 3; c++) {
                jx += K(2, b, c) - K(0, b, c);
                jy += K(b, 2, c) - K(b, 0, c);
                jz += K(b, c, 2) - K(b, c, 0);
            }
        const int sol = sb[i] != 0;
        const double rs = sol ? 1.0 : rho;     /* solidos: evita divisao por valores sem sentido */
        double ux = (jx + 0.5 * Fx) / rs, uy = (jy + 0.5 * Fy) / rs, uz = (jz + 0.5 * Fz) / rs;
        ux = sol ? 0.0 : ux; uy = sol ? 0.0 : uy; uz = sol ? 0.0 : uz;
        ro[i] = sol ? 1.0 : rho; xo[i] = ux; yo[i] = uy; zo[i] = uz;
        fwd_axis(k, uz, 2);
        fwd_axis(k, uy, 1);
        fwd_axis(k, ux, 0);
        const double ir = 1.0 / rs;
        for (int d = 0; d < 27; d++) k[d] *= ir;
        if (model == 0) {
            double k1x = -K(1, 0, 0), k1y = -K(0, 1, 0), k1z = -K(0, 0, 1);
            for (int d = 0; d < 27; d++) k[d] += w1 * (Keq[d] - k[d]);
            K(1, 0, 0) = k1x; K(0, 1, 0) = k1y; K(0, 0, 1) = k1z;
        } else {
            cumulant_relax(k, ux, uy, uz, prm);
        }
        for (int d = 0; d < 27; d++) k[d] *= rs;
        bwd_axis(k, ux, 0);
        bwd_axis(k, uy, 1);
        bwd_axis(k, uz, 2);
        for (int d = 0; d < 27; d++) buf[d][i] = sol ? Wd[d] : k[d];
    }
}

/*
 * Colisao (em todos os nos) + propagacao periodica.
 *   f      : populacoes no instante t (nao sao alteradas)
 *   fout   : recebe as populacoes propagadas (instante t+1).  A pos-colisao do no n na
 *            direcao d fica em fout[d][n + c_d]; as paredes do lbm3d.py a leem de la.
 *   solid  : 1 = no solido (pos-colisao = pesos w, como no lbm3d.py)
 *   force  : forca de corpo (Fx, Fy, Fz)
 *   rho, ux, uy, uz : saida (macroscopicas pre-colisao; solidos: rho = 1, u = 0)
 * Os nos sao processados em blocos (BLK nos de um plano x): leitura e escrita direcao
 * por direcao, em sequencia na memoria (o calculo e limitado pela memoria, nao pela CPU).
 */
void collide_stream(const double *f, double *fout, const uint8_t *solid, int nx, int ny, int nz,
                    const double *force, const double *prm,
                    double *rho_o, double *ux_o, double *uy_o, double *uz_o)
{
    const long N = (long)nx * ny * nz;
    const long P = (long)ny * nz;
    const int model = (int)prm[0];
    const double w1 = prm[1];
    double Wd[27], Keq[27];
    for (int a = 0; a < 3; a++)
        for (int b = 0; b < 3; b++)
            for (int c = 0; c < 3; c++) {
                Wd[a * 9 + b * 3 + c] = W1[a] * W1[b] * W1[c];
                Keq[a * 9 + b * 3 + c] = E1[a] * E1[b] * E1[c];
            }
    const double Fx = force[0], Fy = force[1], Fz = force[2];
    static double buf[27][BLKPAD];
    static long dest[BLKPAD];

    for (int x = 0; x < nx; x++) {
        const long base = (long)x * P;
        for (long p0 = 0; p0 < P; p0 += BLK) {
            const int nb = (int)((P - p0) < BLK ? (P - p0) : BLK);
            for (int d = 0; d < 27; d++) {
                const double *src = f + d * N + base + p0;
                for (int i = 0; i < nb; i++) buf[d][i] = src[i];
            }
            const uint8_t *sb = solid + base + p0;
            double *ro = rho_o + base + p0, *xo = ux_o + base + p0;
            double *yo = uy_o + base + p0, *zo = uz_o + base + p0;
            if (model == 2) bloco(buf, sb, nb, Fx, Fy, Fz, prm, Wd, Keq, ro, xo, yo, zo, 2);
            else if (model == 1) bloco(buf, sb, nb, Fx, Fy, Fz, prm, Wd, Keq, ro, xo, yo, zo, 1);
            else bloco(buf, sb, nb, Fx, Fy, Fz, prm, Wd, Keq, ro, xo, yo, zo, 0);
            for (int a = 0; a < 3; a++) {
                int xn = x + a - 1; xn = xn < 0 ? nx - 1 : (xn >= nx ? 0 : xn);
                for (int b = 0; b < 3; b++)
                    for (int c = 0; c < 3; c++) {
                        const int d = a * 9 + b * 3 + c;
                        double *dst = fout + d * N + (long)xn * P;
                        int y = (int)(p0 / nz), z = (int)(p0 % nz);
                        for (int i = 0; i < nb; i++) {
                            int yn = y + b - 1; yn = yn < 0 ? ny - 1 : (yn >= ny ? 0 : yn);
                            int zn = z + c - 1; zn = zn < 0 ? nz - 1 : (zn >= nz ? 0 : zn);
                            dest[i] = (long)yn * nz + zn;
                            if (++z == nz) { z = 0; y++; }
                        }
                        for (int i = 0; i < nb; i++) dst[dest[i]] = buf[d][i];
                    }
            }
        }
    }
}

/*
 * Passo termico: TRT D3Q7 para a energia E = C T (C = capacidade termica volumetrica do
 * no, relativa a do fluido), com equilibrio
 *   g_0^eq = E - T + w0 T,   g_i^eq = w_i [T + (T - Tref) c_i.u / cs2],   w0 = 1/4, w_i = 1/8, cs2 = 1/4
 * de modo que  dE/dt + div((T - Tref) u) = div(k grad T),  k = cs2 (tau- - 1/2)  (por no).
 * Como div(u) ~ 0, equivale a u.grad(T); Tref escolhe onde fica o erro da divergencia
 * discreta de u (grande na interface solido-fluido): ele e proporcional a T - Tref.
 * Ordem das direcoes: 0 repouso, 1 +x, 2 -x, 3 +y, 4 -y, 5 +z, 6 -z.
 *   g      : entra no instante t, sai com as pos-colisao
 *   gout   : recebe as propagadas
 *   ativo  : 1 = no do dominio termico; 0 = fora (pos-colisao = 0)
 *   C, wm, wp : capacidade relativa, 1/tau- e 1/tau+ por no
 *   ux, uy, uz : velocidade do escoamento (0 nos solidos)
 *   T      : saida, temperatura no instante t (0 fora do dominio)
 */
void thermal_collide_stream(double *g, double *gout, const uint8_t *ativo, int nx, int ny, int nz,
                            const double *C, const double *wm, const double *wp,
                            const double *ux, const double *uy, const double *uz, double *T_o,
                            double Tref)
{
    const long N = (long)nx * ny * nz;
    const double w0 = 0.25, wi = 0.125;
    for (int x = 0; x < nx; x++) {
        int xp = x + 1 >= nx ? 0 : x + 1, xm = x - 1 < 0 ? nx - 1 : x - 1;
        for (int y = 0; y < ny; y++) {
            int yp = y + 1 >= ny ? 0 : y + 1, ym = y - 1 < 0 ? ny - 1 : y - 1;
            for (int z = 0; z < nz; z++) {
                int zp = z + 1 >= nz ? 0 : z + 1, zm = z - 1 < 0 ? nz - 1 : z - 1;
                const long n = ((long)x * ny + y) * nz + z;
                double q[7];
                if (!ativo[n]) {
                    for (int i = 0; i < 7; i++) q[i] = 0.0;
                    T_o[n] = 0.0;
                } else {
                    for (int i = 0; i < 7; i++) q[i] = g[i * N + n];
                    double E = q[0] + q[1] + q[2] + q[3] + q[4] + q[5] + q[6];
                    double T = E / C[n];
                    T_o[n] = T;
                    double om = wm[n], op = wp[n];
                    double ueq[3] = {ux[n], uy[n], uz[n]};
                    double ep = wi * T;
                    q[0] -= op * (q[0] - (E - T + w0 * T));
                    for (int ax = 0; ax < 3; ax++) {
                        double gp = q[1 + 2 * ax], gm = q[2 + 2 * ax];
                        double sp = 0.5 * (gp + gm), am = 0.5 * (gp - gm);
                        double em = wi * (T - Tref) * 4.0 * ueq[ax];
                        double ds = op * (sp - ep), da = om * (am - em);
                        q[1 + 2 * ax] = gp - ds - da;
                        q[2 + 2 * ax] = gm - ds + da;
                    }
                }
                for (int i = 0; i < 7; i++) g[i * N + n] = q[i];
                gout[0 * N + n] = q[0];
                gout[1 * N + ((long)xp * ny + y) * nz + z] = q[1];
                gout[2 * N + ((long)xm * ny + y) * nz + z] = q[2];
                gout[3 * N + ((long)x * ny + yp) * nz + z] = q[3];
                gout[4 * N + ((long)x * ny + ym) * nz + z] = q[4];
                gout[5 * N + ((long)x * ny + y) * nz + zp] = q[5];
                gout[6 * N + ((long)x * ny + y) * nz + zm] = q[6];
            }
        }
    }
}
