#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Implementação do cálculo rápido do polinômio
void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx) {
    *px = p.p[p.grau];
    *dpx = 0.0;
    for (int i = p.grau - 1; i >= 0; i--) {
        *dpx = (*dpx) * x + (*px);
        *px = (*px) * x + p.p[i];
    }
}

// Implementação do cálculo lento do polinômio
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

// Implementação do método de Newton-Raphson
real_t newtonRaphson(Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz, CalcPolinomioFunc calcPolinomio) {
    real_t xk = x0, xk1, fx, dfx, erro;
    *it = 0;
    do {
        calcPolinomio(p, xk, &fx, &dfx);
        if (dfx == 0.0) break;
        xk1 = xk - fx / dfx;
        
        if (criterioParada == 1) {
            erro = fabs(xk1 - xk);
        } else if (criterioParada == 2) {
            erro = fabs(fx);
        } else if (criterioParada == 3){
            erro = fabs(xk1 - xk) / DBL_EPSILON;
        }
        
        xk = xk1;
        (*it)++;
    } while (erro > EPS && *it < MAXIT);
    *raiz = xk;
    return erro;
}

// Implementação do método da bissecção
real_t bisseccao(Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, CalcPolinomioFunc calcPolinomio) {
    real_t fa, fb, fm, m, erro, temp;
    *it = 0;
    
    calcPolinomio(p, a, &fa, &temp);
    calcPolinomio(p, b, &fb, &temp);
    if (fa * fb > 0) return -1;
    
    do {
        m = (a + b) / 2.0;
        calcPolinomio(p, m, &fm, &temp);
        
        
        if (fa * fm < 0) {
            b = m;
            fb = fm;
        } else {
            a = m;
            fa = fm;
        }
        
        if (criterioParada == 1) {
            erro = fabs(b - a);
        } else if (criterioParada == 2) {
            erro = fabs(fm);
        } else if (criterioParada == 3) {
            erro = fabs(b - a) / DBL_EPSILON;
        }
        
        (*it)++;
    } while (erro > EPS && *it < MAXIT);
    
    *raiz = m;
    return erro;
}