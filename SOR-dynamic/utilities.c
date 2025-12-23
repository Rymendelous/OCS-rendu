
#include "utilities.h"
#include "blas2_static.h"
#include "iterative.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>

double u_exact(double x, double y) {
  double res = x * x * x * (x - 1.) * (x - 1.) + y * y * y * (y - 1.);
  return res;
}

double b_rhs(double x, double y) {

    double res = -12. * (x * x - x + 1. / 6. + y * y - 0.5 * y);
  return res;
}

void init_a(double a[Nx][Ny][5], double hx, double hy) {
  size_t i = 0, j = 0, k = 0;
  // boundaries treatment
  for (k = 0; k < 5; k++) {
    // x = 0 or x= 1
    for (j = 0; j < Ny; j++) {
      a[0][j][k] = 0;
      a[Nx - 1][j][k] = 0;
    }

    // y = 0 or y = 1
    for (i = 0; i < Nx; i++) {
      a[i][0][k] = 0;
      a[i][Ny - 1][k] = 0;
    }
  }

  // interior points
  for (i = 1; i <= Nx - 2; i++) {
    for (j = 1; j <= Ny - 2; j++) {
      a[i][j][0] = -1.0 / (hx * hx);
      a[i][j][1] = -1.0 / (hy * hy);
      a[i][j][2] = 2 * (1. / (hx * hx) + 1. / (hy * hy));
      a[i][j][3] = -1.0 / (hx * hx);
      a[i][j][4] = -1.0 / (hy * hy);
    }
  }
}

void init_b(double b[Nx][Ny], double hx, double hy) {

  for (size_t i = 1; i < Nx-1; i++) {
    for (size_t j = 1; j < Ny-1; j++) {
      b[i][j] = b_rhs(i * hx, j * hy);
    }
  }
  for (int i = 0; i < Nx; i++){
    b[i][0] = 0;
    b[i][Ny-1] = 0;
  }

  for (int j = 0; j < Ny; j++){
    b[0][j] = 0;
    b[Nx-1][j] = 0;
  }

}

void set_u(double u[Nx][Ny], double value) {
  for (size_t i = 0; i < Nx; i++) {
    for (size_t j = 0; j < Ny; j++) {
      u[i][j] = value;
    }
  }
}

double l2_nrm_u(double u[Nx][Ny]) {
  double res = 0.;
  for (size_t i = 0; i < Nx; i++) {
    for (size_t j = 0; j < Ny; j++) {
      res += (u[i][j] * u[i][j]);
    }
  }
  return sqrt(res);
}

void copy_u(double src[Nx][Ny], double dst[Nx][Ny]) {

  for (size_t i = 0; i < Nx; i++) {
    for (size_t j = 0; j < Ny; j++) {
      dst[i][j] = src[i][j];
    }
  }
}

// void print_a(double a[Nx][Ny][5]){
//     printf("\n\n");
 
//         printf("\n");
//         for(size_t i = 0; i < Nx; i++){
//             for(size_t j = 0; j < Ny; j++){
//                   for(size_t k =0; k < 5; k++){
//                     printf("%lf   ", a[i][j][k]);
//                   }
//                   printf("\t");
//             }
//             printf("\n"); 
//         }

//     printf("\n\n");
// }

void print_a(double a[Nx][Ny][5]) {
    printf("\n\n");

    for (size_t i = 0; i < Nx; i++) {
        printf("i = %zu:\n", i);

        for (size_t j = 0; j < Ny; j++) {
            printf("  j = %zu :  ", j);

            for (size_t k = 0; k < 5; k++) {
                printf("%10.6f  ", a[i][j][k]);
            }

            printf("\n");
        }

        printf("\n");
    }

    printf("\n\n");
}

void init_a_dyn(double*** a, double hx, double hy) {
    for (size_t i = 0; i < Nx; i++) {
        for (size_t j = 0; j < Ny; j++) {
            for (size_t k = 0; k < 5; k++) {
                if (i == 0 || i == Nx-1 || j == 0 || j == Ny-1) {
                    a[i][j][k] = 0.;
                } else {
                    a[i][j][0] = -1.0 / (hx*hx);
                    a[i][j][1] = -1.0 / (hy*hy);
                    a[i][j][2] = 2*(1/(hx*hx) + 1/(hy*hy));
                    a[i][j][3] = -1.0 / (hx*hx);
                    a[i][j][4] = -1.0 / (hy*hy);
                }
            }
        }
    }
}


void set_u_dyn(double** u, double value) {
    for (size_t i = 0; i < Nx; i++) {
        for (size_t j = 0; j < Ny; j++) {
            u[i][j] = value;
        }
    }
}

void init_b_dyn(double** b, double hx, double hy) {
    // allouer le tableau 2D si ce n'est pas déjà fait
    for (size_t i = 0; i < Nx; i++) {
        b[i] = malloc(Ny * sizeof(double));
    }

    // points intérieurs
    for (size_t i = 1; i < Nx-1; i++) {
        for (size_t j = 1; j < Ny-1; j++) {
            b[i][j] = b_rhs(i * hx, j * hy);
        }
    }

    // conditions aux bords
    for (size_t i = 0; i < Nx; i++) {
        b[i][0] = 0.;
        b[i][Ny-1] = 0.;
    }

    for (size_t j = 0; j < Ny; j++) {
        b[0][j] = 0.;
        b[Nx-1][j] = 0.;
    }
}