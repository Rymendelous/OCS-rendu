#include "iterative.h"
#include "utilities.h"
#include "gs_dynamic.h"
#include "blas1_static.h"
#include "blas2_static.h"
#include <stdlib.h>  // pour malloc et free

int main (int argc, char** argv){

    double*** a_dynamic; 
    double** b_dynamic;
    double** u_old_dynamic;
    double* residual_stories;

    double tol = 1e-10;
    double hx = 1./(Nx-1), hy = 1./(Ny-1);

    a_dynamic = malloc(Nx * sizeof(double**));
    for (size_t i = 0; i < Nx; i++) {
        a_dynamic[i] = malloc(Ny * sizeof(double*));
        for (size_t j = 0; j < Ny; j++) {
            a_dynamic[i][j] = malloc(5 * sizeof(double));
        }
    }

    b_dynamic = malloc(Nx * sizeof(double*));
    for(size_t i = 0;i<Nx;i++){
        b_dynamic[i]=malloc(Ny*sizeof(double));
    }

    u_old_dynamic = malloc(Nx * sizeof(double*));
    for (size_t i = 0; i < Nx; i++) {
        u_old_dynamic[i] = malloc(Ny * sizeof(double));
    }

   residual_stories = malloc(Niter_max * sizeof(double));


    init_a_dyn(a_dynamic, hx, hy);
    init_b_dyn(b_dynamic, hx, hy);
    set_u_dyn(u_old_dynamic, 0.);
    dscal(residual_stories, 0., Niter_max);

   
    gs_dynamic(a_dynamic, b_dynamic, u_old_dynamic, residual_stories, hx, hy, tol);

    
    for (size_t i=0; i < Nx; i++) {
        for (size_t j=0; j < Ny; j++) {
            free(a_dynamic[i][j]);
        }
        free(a_dynamic[i]);
    }

    for(size_t i=0;i<Nx;i++){
        free(b_dynamic[i]);
    }

    for(size_t i=0;i<Nx;i++){
        free(u_old_dynamic[i]);
    }

    free(a_dynamic);
    free(b_dynamic);
    free(u_old_dynamic);
    free(residual_stories);

    return 0;
}
