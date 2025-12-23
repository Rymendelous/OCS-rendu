
#include "iterative.h"
#include "jacobi_static.h"
#include "utilities.h"
#include "blas1_static.h"
#include "blas2_static.h"



int main (int argc, char** argv){

    /*
    static cases
    
    */

    // declarations of variables
    double a_static[Nx][Ny][5]; 
    double b_static[Nx][Nx];
    double u_old_static[Nx][Nx];
    double u_new_static[Nx][Ny];
    double residual_stories[Niter_max];

    double tol = 1e-10;
    double hx = 1./(Nx-1), hy = 1./(Ny -1);

    // initialize a and b
    init_a(a_static, hx, hy);
    // print_a(a_static);

    init_b(b_static, hx, hy);
    // print2D(b_static);

    // initialize u_old and u_new by 0
    set_u(u_old_static, 0.);
    // print2D(u_old_static);
    set_u(u_new_static, 0.);
    // print2D(u_new_static);


    // initialize residual_stories
    dscal(residual_stories, 0., Niter_max);
    // print1D(residual_stories, Niter_max);

    // perform Jacobi iterations
    jacobi_static(a_static, b_static, u_old_static, u_new_static, residual_stories, hx, hy, tol);
    // print1D(residual_stories, Niter_max);
    print2D(u_new_static);
    


    return 0;
}