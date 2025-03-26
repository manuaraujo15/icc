#include <stdio.h>
#include <math.h>

#include "utils.h"
#include "ZeroFuncao.h"

// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz){
real_t px, dpx;
    real_t erro = 1.0;
    *it = 0;
    
    while (erro > criterioParada) {
        calcPolinomio_rapido(p, x0, &px, &dpx); // Calcula valor e derivada
        if (dpx == 0) {
            // Evitar divisão por zero
            return -1;
        }
        real_t x1 = x0 - px / dpx; // Fórmula de Newton-Raphson
        erro = fabs(x1 - x0); // Erro absoluto
        x0 = x1;
        (*it)++;
    }
    *raiz = x0;
    return erro;
}


// Retorna valor do erro quando método finalizou. Este valor depende de tipoErro
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz){
    real_t fa, fb, c, fc;
    *it = 0;
    calcPolinomio_rapido(p, a, &fa, &fb); // Calcula f(a) e f(b)
    
    if (fa * fb > 0) {
        return -1; // Não há garantia de raiz no intervalo [a, b]
    }
    
    while ((b - a) > criterioParada) {
        c = (a + b) / 2;
        calcPolinomio_rapido(p, c, &fc, &fc); // Calcula f(c)
        
        if (fc == 0.0) {
            *raiz = c;
            return 0.0;
        }
        
        if (fa * fc < 0) {
            b = c;
        } else {
            a = c;
            fa = fc;
        }
        (*it)++;
    }
    
    *raiz = (a + b) / 2;
    return fabs(b - a);
}


void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx){
    *px = 0.0;
    *dpx = 0.0;
    real_t potencia = 1.0;
    
    for (int i = 0; i <= p.grau; i++) {
        *px += p.p[i] * potencia;
        if (i > 0) {
            *dpx += p.p[i] * i * potencia / x;
        }
        potencia *= x;
    }
}


void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx){
    *px = 0.0;
    *dpx = 0.0;
    
    for (int i = 0; i <= p.grau; i++) {
        *px += p.p[i] * pow(x, i);
        if (i > 0) {
            *dpx += p.p[i] * i * pow(x, i - 1);
        }
    }
}
