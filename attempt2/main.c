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
//    PrintMat(G, n, n);
    //printf("\n");
    U = G;
    //PrintMat(U, n, n);

    GET_TIME(start);
    if (n == 1)
        x[0] = U[0][1] / U[0][0];
    else {
//# pragma omp parallel num_threads(thread_count) \
//          shared(U, n, row) private(i, j, temp, k)
//        {
//            omp_set_num_threads(1);
//            printf("%d\n", omp_get_num_threads( ));

//# pragma omp sections nowait
//            {
//# pragma omp section
//printf("%d\n", thread_count);
                for (k = 0; k < n - 1; ++k) {
                    temp = 0;
                    j = 0;
//                    omp_set_num_threads(thread_count);
//# pragma omp for schedule(static)

//                    # pragma omp section
                        for (i = k; i < n; ++i) {
                            // if the value is smaller, change temp to that value.
//# pragma omp ordered

                            if (temp < U[row[i]][k] * U[row[i]][k]) {
                                temp = U[row[i]][k] * U[row[i]][k];
                                j = i;
                            }
                        }
//#pragma omp barrier
//#pragma omp master
//                    {printf("%d\n", omp_get_num_threads( ));}
//                    omp_set_num_threads(1);
//                    # pragma omp section
                        if (j != k) {
                            i = row[j];
                            row[j] = row[k];
                            row[k] = i;
                            /* row replacement*/
                        }
                        /*calculating*/
//# pragma omp for schedule(static)
//                    # pragma omp section
//# pragma omp barrier
//# pragma omp for schedule(dynamic,(n/thread_count)) nowait
                # pragma omp parallel for num_threads(thread_count) default(none) shared(U, row, k, n) private(i, j, temp)
                        for (i = k + 1; i < n; ++i) {
                            temp = U[row[i]][k] / U[row[k]][k];
//                    # pragma omp for
                            for (j = k; j < n + 1; ++j)
//# pragma omp critical
                                    U[row[i]][j] -= U[row[k]][j] * temp;
                        }
//# pragma omp barrier

//                }

                }

                /*Jordan elimination*/

//# pragma omp section
//# pragma omp parallel for num_threads(thread_count) collapse(2)
//# pragma omp parallel for num_threads(thread_count) shared(U, row, k, n) private(i, temp)
//        {
//#pragma omp for collapse(2)
//            omp_set_num_threads(1);
        for (k = n - 1; k > 0; --k) {
//            omp_set_num_threads(thread_count);
//# pragma omp for
            for (i = k - 1; i >= 0; --i) {
                        temp = U[row[i]][k] / U[row[k]][k];
//                    # pragma omp critical
                        U[row[i]][k] -= temp * U[row[k]][k];
//                    # pragma omp critical
                        U[row[i]][n] -= temp * U[row[k]][n];
                    }
                }
//# pragma omp section
//#pragma omp for
                for (k = 0; k < n; ++k) {
//            printf("6 \n");
                    x[k] = U[row[k]][n] / U[row[k]][k];
                }
//            }
//        }
    }

    D = U;
    GET_TIME(end);
    time = end-start;
    Lab3SaveOutput(x, n, time);
    PrintMat(D, n, n);
//    printf("\n");
//    PrintVec(x, n);

    return 0;
}
