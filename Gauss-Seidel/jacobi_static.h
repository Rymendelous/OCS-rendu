#pragma once
#include <stddef.h>
#include<stdio.h>
#include<math.h>
#include "iterative.h"

void jacobi_static(double a[Nx][Ny][5], double b[Nx][Ny], double u_old[Nx][Ny], double u_new[Nx][Ny], double residual [Niter_max], double hx, double hy, double tol);

// void jacobi_dynamic(double*** a, double** u_old, double** u_new, double* residual , double hx, double hy, double tol);