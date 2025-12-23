#pragma once
#include <stddef.h>
#include<stdio.h>
#include<math.h>
#include "iterative.h"

void jacobi_dynamic(double*** a, double** u_old, double** u_new, double* residual , double hx, double hy, double tol);