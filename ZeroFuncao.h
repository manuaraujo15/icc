#ifndef __ZEROFUNCAO_H__
#define __ZEROFUNCAO_H__

#include <float.h>

// Aproximação aceitável como valor zero
#define ZERO DBL_EPSILON

// Parâmetros para teste de convergência
#define MAXIT 500
#define EPS 1.0e-6
#define ULPS 2

typedef struct {
  real_t *p;
  int grau;
} Polinomio;


//COLOCA OQ CADA VARIAVEL FAZ
/*
  p é um vetor de tamanho grau+1 
  os parâmetros a, b, x0 são os chutes iniciais para os métodos,
  it é o número de iterações do método que foram necessários para chegar à raiz, cujo valor deve ser devolvido pelo parâmetro raiz. 
*/

// Métodos
// Retornam valor do erro quando método finalizou. Este valor depende de tipoErro

real_t newtonRaphson (Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz);
real_t bisseccao (Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz);

// Cálculo de Polinômios
//cálculo do valor px de um polinômio p e sua primeira derivada dpx em um ponto x.
void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx );
void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx );

#endif // __ZEROFUNCAO_H__

