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
    int thread_count = strtol(argv[1], NULL, 10);
    int k, i, j, tid;
    double temp;
    int *row;


    Lab3LoadInput(&G, &n);
    b = CreateVec(n);
    x = CreateVec(n);
    row = malloc(n * sizeof(int));

    #pragma omp parallel for
    for (i = 0; i < n; ++i)
        row[i] = i;

    U = CreateMat(n, n);
    U = G;

    GET_TIME(start);
    if (n == 1)
        x[0] = U[0][1] / U[0][0];
    else {
        for (k = 0; k < n - 1; ++k) {
            temp = 0;
            j = 0;

                for (i = k; i < n; ++i) {
                    // if the value is smaller, change temp to that value.
                    if (temp < U[row[i]][k] * U[row[i]][k]) {
                        temp = U[row[i]][k] * U[row[i]][k];
                        j = i;
                    }
                }
                if (j != k) {
                    i = row[j];
                    row[j] = row[k];
                    row[k] = i;
                    /* row replacement*/
                }
                /*calculating*/
            #pragma omp parallel for num_threads(thread_count) default(none) shared(U, row, k, n) private(i, j, temp)
                for (i = k + 1; i < n; ++i) {
                    temp = U[row[i]][k] / U[row[k]][k];
                    for (j = k; j < n + 1; ++j)
                            U[row[i]][j] -= U[row[k]][j] * temp;
                }
        }

                /*Jordan elimination*/
        for (k = n - 1; k > 0; --k) {
            for (i = k - 1; i >= 0; --i) {
                        temp = U[row[i]][k] / U[row[k]][k];
                        U[row[i]][k] -= temp * U[row[k]][k];
                        U[row[i]][n] -= temp * U[row[k]][n];
            }
        }
        for (k = 0; k < n; ++k) {
            x[k] = U[row[k]][n] / U[row[k]][k];
        }
    }

    D = U;
    GET_TIME(end);
    time = end-start;
    Lab3SaveOutput(x, n, time);
    PrintMat(D, n, n);

    return 0;
}
