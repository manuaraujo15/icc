#include <stdio.h>
#include <math.h>
#include <float.h>
#include <fenv.h>
#include "utils.h"
#include "ZeroFuncao.h"
//funcao polinomio
int main() {
  fesetround(FE_DOWNWARD);
  int grau, it;
  real_t a, b, raiz, erro, tempo;
  
  scanf("%d", &grau);
  Polinomio p;
  p.grau = grau;
  p.p = (real_t *)malloc((grau + 1) * sizeof(real_t));
  for (int i = grau; i >= 0; i--) {
      scanf("%lf", &p.p[i]);
  }
  scanf("%lf %lf", &a, &b);
  
  printf("RAPIDO\n");
  for (int criterio = 1; criterio <= 3; criterio++) {
      tempo = timestamp();
      erro = bisseccao(p, a, b, criterio, &it, &raiz);
      tempo = timestamp() - tempo;
      printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
  }
  for (int criterio = 1; criterio <= 3; criterio++) {
      tempo = timestamp();
      erro = newtonRaphson(p, (a + b) / 2, criterio, &it, &raiz);
      tempo = timestamp() - tempo;
      printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
  }
  
  printf("LENTO\n");
  for (int criterio = 1; criterio <= 3; criterio++) {
      tempo = timestamp();
      erro = bisseccao(p, a, b, criterio, &it, &raiz);
      tempo = timestamp() - tempo;
      printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
  }
  for (int criterio = 1; criterio <= 3; criterio++) {
      tempo = timestamp();
      erro = newtonRaphson(p, (a + b) / 2, criterio, &it, &raiz);
      tempo = timestamp() - tempo;
      printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
  }
  
  free(p.p);
  return 0;
}


