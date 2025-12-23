
#include "jacobi_static.h"
#include "iterative.h"
#include "utilities.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>

void jacobi_static(double a[Nx][Ny][5], double b[Nx][Ny], double u_old[Nx][Ny], double u_new[Nx][Ny], double residual[Niter_max], double hx, double hy, double tol){

    size_t it = 0; 
    double cur_res = 1;
    for (it = 0; it < Niter_max; it++){
        printf("it : %ld \t res = %lf \n", it, cur_res);
       

       double   diff_acc = 0.;

        for(size_t i = 0; i < Nx; i++){
            for(size_t j =0; j< Ny; j++){

                // boundaries treatment
                if(i == 0  || i == Nx-1 || j==0 || j== Ny-1)
                    u_new[i][j] = b[i][j];
                else{
                    u_new[i][j] = (b[i][j] - a[i][j][0]* u_old[i-1][j] -a[i][j][1]*u_old[i][j-1] - a[i][j][3]*u_old[i][j+1] - a[i][j][4]*u_old[i+1][j])/a[i][j][2];
                }
                diff_acc += (u_new[i][j] -  u_old[i][j])*((u_new[i][j] -  u_old[i][j]));
            }
        }
        cur_res = sqrt(diff_acc);
        residual[it] = cur_res;

        copy_u(u_new, u_old);


        if(cur_res < tol && it >= 2){
            break;
        }

    }

    if(it < Niter_max)
        printf("[Jacobi] converged after %ld iterations !\n", it);
    else
        printf("[Jacobi] maximum iteration reached without converged !\n");
}


/*

void jacobi_dynamic(double*** a, double** u_old, double** u_new, double* residual , double hx, double hy, double tol){
    size_t it; 
    double cur_res = 1;
    for (it = 0; it < Niter_max; it++){
       

       double   diff_acc = 0.;

        for(size_t i = 0; i < Nx; i++){
            for(size_t j =0; j< Ny; j++){

                double x_bord = (i==0) ? 0 : 1;
                double y_bord = (j==0) ? 0 : 1;

                // boundaries treatment
                if(i == 0  || i == Nx-1 || j==0 || j== Ny-1)
                    u_new[i][j] = b_rhs(i*hx, j*hy, x_bord, y_bord);
                else{
                    u_new[i][j] = (b_rhs(i*hx, j*hy, x_bord, y_bord) - a[i][j][0]* u_old[i-1][j] -a[i][j][1]*u_old[i][j-1] - a[i][j][3]*u_old[i][j+1] - a[i][j][4]*u_old[i+1][j])/a[i][j][2];
                }

                diff_acc += (u_new[i][j] -  u_old[i][j]);
            }
        }
        cur_res = sqrt(diff_acc);
        residual[it] = cur_res;
        it++;

        copy_u_dyn(u_new, u_old);


        if(fabs(cur_res - tol) < 0){
            break;
        }

    }

    if(it < Niter_max)
        printf("[Jacobi] converged after %ld iterations !\n", it);
    else
        printf("[Jacobi] maximum iteration reached without converged !\n");
}

*/