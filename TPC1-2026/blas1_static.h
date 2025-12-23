
#pragma once
#include <math.h>
#include <stddef.h>
#include <stdio.h>

//%%%%%%%%%%%%%%%%%%%%%%%%%%%% declaration [functions] %%%%%%%%%%%%%%%%

void ddot(double v1[], double v2[], double *res, size_t size);

void daxpy(double alpha, double x[], double y[], size_t size);

void print1D(double v[], size_t size);

void dscal(double v[], double alpha, size_t size);

void copy_v(double src[], double dst[], size_t size);

void l2norm(double v[], double *res, size_t size);

void l1norm(double v[], double *res, size_t size);

