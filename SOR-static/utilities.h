#pragma once
#include<stdio.h>
#include<math.h>
#include "iterative.h"


double u_exact(double x, double y);
double b_rhs(double x, double y);
void  init_a(double a[Nx][Ny][5], double hx, double hy);
void init_b(double b[Nx][Ny], double hx, double hy);
void set_u(double u[Nx][Ny], double value);
double l2_nrm_u(double u[Nx][Ny]);
void copy_u(double src[Nx][Ny], double dst[Nx][Ny]);
void print_a(double a[Nx][Ny][5]);

void  init_a_dyn(double*** a, double hx, double hy);
void set_u_dyn(double** u, double value);
double l2_nrm_u_dyn(double** u);
void copy_u_dyn(double** src, double** dst);




