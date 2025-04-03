#include <stdio.h>
#include <math.h>
#include <float.h>

#include "utils.h"
#include "ZeroFuncao.h"
#include <fenv.h>
//funcao polinomio
int main() {
    fesetround(FE_DOWNWARD);  // Define arredondamento para baixo
    
    int grau;
    scanf("%d", &grau);
    
    Polinomio p;
    p.grau = grau;
    p.p = (real_t*) malloc((grau + 1) * sizeof(real_t));
    if (!p.p) {
        printf("Erro ao alocar memória!\n");
        return -1;
    }                                                                                                                                                           
    for (int i = 0; i <= grau; i++) {
        scanf("%lf", &p.p[i]);
    }
    
    real_t a, b;
    scanf("%lf %lf", &a, &b);
    
    printf("RAPIDO\n");
    for (int criterio = 1; criterio <= 3; criterio++) {
        int it;
        real_t raiz, erro;
        double tempo;
        tempo = timestamp();
        erro = bisseccao(p, a, b, criterio, &it, &raiz, calcPolinomio_rapido);
        tempo = timestamp() - tempo;
        printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
    }
    for (int criterio = 1; criterio <= 3; criterio++) {
        int it;
        real_t raiz, erro; 
        double tempo;

        tempo = timestamp();
        erro = newtonRaphson(p, (a + b) / 2.0, criterio, &it, &raiz, calcPolinomio_rapido);
        tempo = timestamp() - tempo;
        printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
    }
    
    printf("LENTO\n");
    for (int criterio = 1; criterio <= 3; criterio++) {
        int it;
        real_t raiz, erro;
        double tempo;

        tempo = timestamp();
        erro = bisseccao(p, a, b, criterio, &it, &raiz, calcPolinomio_lento);
        tempo = timestamp() - tempo;
        printf("bissec %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
    }
    for (int criterio = 1; criterio <= 3; criterio++) {
        int it;
        real_t raiz, erro;
        double tempo;

        tempo = timestamp();
        erro = newtonRaphson(p, (a + b) / 2.0, criterio, &it, &raiz, calcPolinomio_lento);
        tempo = timestamp() - tempo;
        printf("newton %.15e %.15e %d %.8e\n", raiz, erro, it, tempo);
    }
    
    free(p.p);
    return 0;
}