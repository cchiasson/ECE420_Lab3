#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab3IO.h"
#include "timer.h"
#include "Lab3IO.c"
#include "Lab3IO.h"
#include <omp.h>

double **G;
double **U;
double **D;
double* b;
int n;
double* x;
double time, start, end;

int main(int argc, char* argv[]) {
    int thread_count = atoi(argv[1]);
    int k, i, j;
    double temp;
    int *row;


    Lab3LoadInput(&G, &n);
    b = CreateVec(n);
    x = CreateVec(n);
    row = malloc(n * sizeof(int));
    for (i = 0; i < n; ++i)
        row[i] = i;

    U = CreateMat(n, n);
//    PrintMat(G, n, n);
    //printf("\n");
    U = G;
    //PrintMat(U, n, n);

    GET_TIME(start);
    if (n == 1)
        x[0] = U[0][1] / U[0][0];
    else {
#       pragma omp parallel num_threads(thread_count) \
          default(none) shared(U, row, n, k) private(i, j, temp)
        for (k = 0; k < n - 1; ++k) {
            temp = 0;
//            printf("1 \n");
#            pragma omp for
            for (i = k, j = 0; i < n; ++i) {
//                printf("2 \n");
                // if the value is smaller, change temp to that value.
                if (temp < U[row[i]][k] * U[row[i]][k]) {
//                    printf("3 \n");
                    temp = U[row[i]][k] * U[row[i]][k];
                    j = i;
                }
            }

            if (j != k) {
//                printf("3.5 \n");
                i = row[j];
                row[j] = row[k];
                row[k] = i;
                /* row replacement*/
            }
            /*calculating*/
#            pragma omp for
//#            pragma omp parallel for num_threads(thread_count) default(none) shared(U, row, k, n) private(i, j, temp)
            for (i = k + 1; i < n; ++i) {
                temp = U[row[i]][k] / U[row[k]][k];
                for (j = k; j < n + 1; ++j)
                    U[row[i]][j] -= U[row[k]][j] * temp;
            }

        }

#       pragma omp parallel num_threads(thread_count) \
          default(none) shared(U, row, n, k) private(i, temp)
        /*Jordan elimination*/
        for (k = n - 1; k > 0; --k) {
#            pragma omp for
            for (i = k - 1; i >= 0; --i) {
                temp = U[row[i]][k] / U[row[k]][k];
                U[row[i]][k] -= temp * U[row[k]][k];
                U[row[i]][n] -= temp * U[row[k]][n];
            }
        }

#       pragma omp for num_threads(thread_count) \
          default(none) shared(U, row, n, k) private()
        for (k = 0; k < n; ++k) {
//            printf("6 \n");
            x[k] = U[row[k]][n] / U[row[k]][k];
        }



    }

    D = U;
    GET_TIME(end);
    time = end-start;
    Lab3SaveOutput(x, n, time);
//    PrintMat(D, n, n);
//    printf("\n");
//    PrintVec(x, n);

    return 0;
}
