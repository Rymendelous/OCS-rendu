#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include "blas1_static.h"

// int main(int argc, char **argv) {

//   int N1 = 5;

//   double x1[] = {1., 2., -1., 0., 1};
//   double x2[5] = {1., 2., -2., 1., 0.};
//   double r_ddot = 0.0;
//   double alpha = 0.5;
//   double l1nrm = 0.;
//   double l2nrm = 0.;

//   ddot(x1, x2, &r_ddot, N1);
//   printf("<x1, x2 > = %lf", r_ddot);

//   daxpy(alpha, x1, x2, N1);
//   print1D(x2, N1);

//   dscal(x1, 2, N1);
//   print1D(x1, N1);

//   l1norm(x2, &l2nrm, N1);
//   printf("\nL2 norm of x2 = %lf", l2nrm);
//   l2norm(x2, &l1nrm, N1);
//   printf("\nL1 norm of x2 = %lf", l1nrm);

//   return 0;
// }

//%%%%%%%%%%%%%%%%%%%%%%%%%%%% definition [functions] %%%%%%%%%%%%%%%%
void ddot(double v1[], double v2[], double *res, size_t size) {

  *res = 0.;
  for (size_t i = 0; i < size; i++) {
    *res = *res + v1[i] * v2[i];
  }
}

void daxpy(double alpha, double x[], double y[], size_t size) {
  for (size_t i = 0; i < size; i++) {
    y[i] = alpha * x[i] + y[i];
  }
}

void print1D(double v[], size_t size) {
  printf("\n");
  for (size_t i = 0; i < size; i++) {
    printf("%lf \t", v[i]);
  }
  printf("\n");
}

void dscal(double v[], double alpha, size_t size) {
  for (size_t i = 0; i < size; i++) {
    v[i] = alpha * v[i];
  }
}

void copy_v(double src[], double dst[], size_t size) {
  for (size_t i = 0; i < size; i++) {
    dst[i] = src[i];
  }
}

void l2norm(double v[], double *res, size_t size) {
  ddot(v, v, res, size);
  *res = sqrt(*res);
}

void l1norm(double v[], double *res, size_t size) {
  *res = 0;
  for (size_t i = 0; i < size; i++) {
    *res += fabs(v[i]);
  }
}