#pragma once
#include <math.h>
#include <stddef.h>
#include <stdio.h>
#include "iterative.h"

// #define N 4
// #define M 4

//%%%%%%%%%%%%%%%%%%%%%%%%%%%% declaration [functions] %%%%%%%%%%%%%%%%
void dgemv(double alpha, double A[Nx][Ny], double x[Ny], double beta, double y[Nx]);

void print2D(double A[Nx][Ny]);