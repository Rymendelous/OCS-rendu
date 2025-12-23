#pragma once
#include <stddef.h>
#include<stdio.h>
#include<math.h>
#include "iterative.h"

void sor_dynamic(double*** a, double** b, double** u_old, double* residual, double hx, double hy, double tol,double omega);