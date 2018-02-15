#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include "Lab3IO.h"
#include "timer.h"
#include "Lab3IO.c"

double **G;
int **U;
int **D;
double* b;
int n;
double* x;
double time;

int Lab3LoadInput();
int Lab3SaveOutput();
double** CreateMat();
int DestroyMat();
int PrintMat();
double* CreateVec();
int PrintVec();
int DestroyVec();

int main(int argc, char* argv[]) {
  int num_threads = atoi(argv[1]);
  int k, i, j, temp;
  int* row;
  

  
  Lab3LoadInput(&G, &n);
  b = CreateVec(n);
  x = CreateVec(n);
  row = malloc(n * sizeof(int));
  for (i = 0; i < n; ++i)
    row[i] = i;

  U = CreateMat(n, n);
  //PrintMat(G, n, n);
  //printf("\n");
  U = G;
  //PrintMat(U, n, n);
	
  if (n == 1)
        x[0] = U[0][1] / U[0][0];
  else{
  for (k=0; k < n-1; ++k) {
      temp = 0;
    printf("1 \n");
    for (i = k, j=0; i < n; ++i) {
      printf("2 \n");
	// if the value is smaller, change temp to that value.
        if (temp < abs(U[row[i]][k])) {
          printf("3 \n");
          temp = abs(U[row[i]][k]);
          j = i;
        }
    }
    /*calculating*/
    for (i = k + 1; i < n; ++i){
        temp = U[row[i]][k] / U[row[k]][k];
        for (j = k; j < n + 1; ++j)
            U[row[i]][j] -= U[row[k]][j] * temp;
    } 

    // if they are not the same row, swap them.
    if (j != k) {
      printf("3.5 \n");
      i = row[j];
      row[j] = row[k];
      row[k] = i;
      /* row replacement*/
    }

  } 

  /**for (k = n; k > 2; k--) {
    //eliminate elements to zero for each column one after another
    //printf("4 \n");
    for (i = 1; i < k-1; i++) {
      // replace rows one row after another 
      // value in row i/ value in row k. 
      printf("5 \n");
      temp = U[row[i]][k] / U[row[k]][k];
      U[row[i]][k] -= temp * U[row[k]][k];
      U[row[i]][n] -= temp * U[row[k]][n];
    }
  }*/

   /*Jordan elimination*/
  for (k = n - 1; k > 0; --k){
    for (i = k - 1; i >= 0; --i ){
        temp = U[row[i]][k] / U[row[k]][k];
        U[row[i]][k] -= temp * U[row[k]][k];
        U[row[i]][n] -= temp * U[row[k]][n];
    } 
  }

  

  for (k=0; k< n; ++k)
    printf("6 \n");
    x[k] = U[row[k]][n] / U[row[k]][k];
  }

  D = U;
  //Lab3SaveOutput(&x, n, time);

  
  return 0;
}
