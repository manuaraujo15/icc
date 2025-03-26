#include <stdio.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "ZeroFuncao.h"
//funcao polinomio
int main (){

  real_t a, b;
  Polinomio pol;
  real_t *px, *dpx;
  scanf("%d", &pol.grau);

  for (int i=pol.grau; i >=0; --i)
    scanf("%lf", &pol.p[i]);

  scanf("%lf %lf", &a, &b); // intervalo onde está uma das raizes.
  calcPolinomio_rapido(pol, a, &px, &dpx); // Calcula f(a) e f(b)
  calcPolinomio_lento(pol, a,&px, &dpx );


  // Restante do programa a partir daqui

  return 0;
}

