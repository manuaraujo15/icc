#ifndef __ZEROFUNCAO_H__
#define __ZEROFUNCAO_H__

#include <float.h>
#include "utils.h"

// Aproximação aceitável como valor zero
#define ZERO DBL_EPSILON

// Parâmetros para teste de convergência

#define ZERO DBL_EPSILON
#define MAXIT 500
#define EPS 1.0e-6
#define ULPS 2
typedef double real_t; // Definição do tipo real_t como um alias para double

// Estrutura que representa um polinômio
typedef struct {
    real_t *p;  // Vetor de coeficientes do polinômio
    int grau;    // Grau do polinômio
} Polinomio;

// Tipo de função para cálculo do polinômio e sua derivada
typedef void (*CalcPolinomioFunc)(Polinomio, real_t, real_t*, real_t*);

// Métodos numéricos para encontrar raízes de polinômios
// Retornam o erro absoluto |f(xk)| da última iteração
// Parâmetros:
// - p: Estrutura contendo os coeficientes e o grau do polinômio
// - x0: Estimativa inicial da raiz (para Newton-Raphson)
// - a, b: Intervalo inicial (para bissecção)
// - criterioParada: Define o critério de parada do método
// - it: Ponteiro para armazenar o número de iterações realizadas
// - raiz: Ponteiro para armazenar a raiz encontrada
// - calcPolinomio: Função para calcular o polinômio e sua derivada
// Retorno: erro absoluto da última iteração
real_t newtonRaphson(Polinomio p, real_t x0, int criterioParada, int *it, real_t *raiz, CalcPolinomioFunc calcPolinomio);
real_t bisseccao(Polinomio p, real_t a, real_t b, int criterioParada, int *it, real_t *raiz, CalcPolinomioFunc calcPolinomio);

// Cálculo de Polinômios
// Calcula o valor do polinômio e sua derivada para um dado x
void calcPolinomio_rapido(Polinomio p, real_t x, real_t *px, real_t *dpx);
void calcPolinomio_lento(Polinomio p, real_t x, real_t *px, real_t *dpx);

#endif // __ZEROFUNCAO_H__
