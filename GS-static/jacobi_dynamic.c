
#include "jacobi_static.h"
#include "iterative.h"
#include "utilities.h"
#include <math.h>
#include <stddef.h>
#include <stdio.h>


void jacobi_dynamic(double*** a, double** b, double** u_old, double** u_new, double* residual , double hx, double hy, double tol){
    size_t it; 
    double cur_res = 1;

    for (it = 0; it < Niter_max; it++){
       

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
        it++;

        copy_u_dyn(u_new, u_old);


        if(cur_res < tol && it >=2){
            break;
        }

    }

    if(it < Niter_max)
        printf("[Jacobi-dynamic] converged after %ld iterations !\n", it);
    else
        printf("[Jacobi-dynamic] maximum iteration reached without converged !\n");
}
