
#include <stddef.h>
#include <stdio.h>
#include "blas2_static.h"


// int main(int argc, char **argv) {
//   double A[N][M] = {{1.0, 2.0, 3.0, 1.0},
//                     {2.0, -5., -0.25, 1.},
//                     {3.0, -1.0, 2.0, 1.0},
//                     {-1.0, 2.0, -3.0, 1.}};

//   print2D(A);

//   double x[] = {2.0, -5., -0.25, 1.};
//   double y[] = {1.0, 2.0, 3.0, 1.0};
//   double alpha = 0.5, beta = -1;

//   dgemv(alpha, A, x, beta, y);

//   // TODO print y

//   return 0;
// }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%% definition [functions] %%%%%%%%%%%%%%%%
void dgemv(double alpha, double A[Nx][Ny], double x[Ny], double beta,
           double y[Nx]) {

  for (size_t i = 0; i < Nx; i++) {
    double ligne_sum = 0.0;
    for (size_t j = 0; j < Ny; j++) {
      ligne_sum += A[i][j] * x[j];
    }
    y[i] = ligne_sum;
  }
}

void print2D(double A[Nx][Ny]) {
  printf("\n");

  for (size_t i = 0; i < Nx; i++) {
    printf("\n");
    for (size_t j = 0; j < Ny; j++) {
      printf("%lf \t", A[i][j]);
    }
    printf("\n");
  }

  printf("\n");
}