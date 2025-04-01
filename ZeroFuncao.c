#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx) {
    *px = p.p[p.grau];
    *dpx = 0.0;
    for (int i = p.grau - 1; i >= 0; i--) {
        *dpx = (*dpx) * x + (*px);
        *px = (*px) * x + p.p[i];
    }
}

void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx) {
    *px = 0.0;
    *dpx = 0.0;
    for (int i = 0; i <= p.grau; i++) {
        *px += p.p[i] * pow(x, i);
        if (i > 0) {
            *dpx += i * p.p[i] * pow(x, i - 1);
        }
    }
}

real_t bisseccao(Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz) {
    real_t fa, fb, fc, c;
    int iter = 0;
    do {
        c = (a + b) / 2;
        calcPolinomio_rapido(p, a, &fa, &fc);
        calcPolinomio_rapido(p, c, &fc, &fb);
        if (fa * fc < 0) {
            b = c;
        } else {
            a = c;
        }
        iter++;
    } while (fabs(b - a) > 1e-6 && iter < 500);
    *raiz = c;
    *it = iter;
    return fabs(fc);
}

real_t newtonRaphson(Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz) {
    real_t xk = x0, fx, dfx;
    int iter = 0;
    do {
        calcPolinomio_rapido(p, xk, &fx, &dfx);
        if (fabs(dfx) < DBL_EPSILON) break;
        xk -= fx / dfx;
        iter++;
    } while (fabs(fx) > DBL_EPSILON && iter < 500);
    *raiz = xk;
    *it = iter;
    return fabs(fx);
}