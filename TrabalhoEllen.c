#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define DADOS 100

double N[DADOS];
double t[DADOS];
int n = 0;

int main(int argc, char *argv[]) {
    if (argc != 2) {
        printf("Uso: %s nome_do_arquivo\n", argv[0]);
        return 1;
    }

    FILE *file = fopen(argv[1], "r");

    if (file == NULL) {
        printf("Erro ao abrir o arquivo.\n");
        return 1;
    }

    while (fscanf(file, "%lf %lf", &t[n], &N[n]) == 2) {
        n++;

        if (n >= DADOS) {
            printf("Numero maximo de pontos de dados excedido.\n");
            break;
        }
    }

    fclose(file);

    double beta0, beta1, beta2;
    double somaN = 0.0, somaT = 0.0, somaNT = 0.0, somaT2 = 0.0, somaT3 = 0.0, somaNT2 = 0.0;

    for (int i = 0; i < n; i++) {
        double ti = t[i];
        double ti2 = ti * ti;

        somaN += N[i];
        somaT += ti;
        somaNT += N[i] * ti;
        somaT2 += ti2;

        somaT3 += ti * ti2;
        somaNT2 += N[i] * ti2;
    }

    double denominador1 = n * somaT2 - somaT * somaT;
    beta0 = (somaT2 * somaN - somaT * somaNT) / denominador1;
    beta1 = (n * somaNT - somaT * somaN) / denominador1;

    double denominador2 = n * somaT3 - somaT2 * somaT;
    beta2 = (n * somaNT2 - somaNT * somaT) / denominador2;

    double r1 = 0.0, r2 = 0.0, r3 = 0.0; //coeficiente de determinação (R²)
    double sqr1 = 0.0, sqr2 = 0.0, sqr3 = 0.0, stq = 0.0;

    for (int i = 0; i < n; i++) {
        double estimado1 = beta0 + beta1 * t[i];
        double estimado2 = beta0 + beta1 * t[i] + beta2 * t[i] * t[i];
        double estimado3 = beta0 * exp(beta1 * t[i]);

        sqr1 += (estimado1 - somaN / n) * (estimado1 - somaN / n);
        sqr2 += (estimado2 - somaN / n) * (estimado2 - somaN / n);
        sqr3 += (estimado3 - somaN / n) * (estimado3 - somaN / n);
        stq += (N[i] - somaN / n) * (N[i] - somaN / n);
    }

    r1 = 1.0 - sqr1 / stq; //somas dos quadrados residuais r1 /soma dos quadrados totais r1
    r2 = 1.0 - sqr2 / stq; //somas dos quadrados residuais r2 /soma dos quadrados totais r2
    r3 = 1.0 - sqr3 / stq; //somas dos quadrados residuais r3 /soma dos quadrados totais r3

    if (r1 > r2 && r1 > r3) {
        printf("O Modelo 1 tem o melhor ajuste com um coeficiente de determinacao de %.4f\n", r1);
    } else if (r2 > r1 && r2 > r3) {
        printf("O Modelo 2 tem o melhor ajuste com um coeficiente de determinacao de %.4f\n", r2);
    } else {
        printf("O Modelo 3 tem o melhor ajuste com um coeficiente de determinacao de %.4f\n", r3);
    }

    return 0;
}
