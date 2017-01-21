#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "slu_ddefs.h"
#include "../superlu_util.h"

#define LOC(i,j,k) ( ((i+nx)%nx) + nx*((j+ny)%ny) + nx*ny*((k+nz)%nz) )

#define alpha (1./2. )
#define beta  (1./20.)
#define a10   ((17./12.)  /2.)
#define b10   ((101./150.)/4.)
#define c10   ((1./100.)  /6.)

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

void create_pade_matrix(double *nzval, int *rowind, int *colptr, long int nnz, int nx, int ny, int nz)
{
    // double *nzval;
    // int *rowind, *colptr;
    // long int nnz;

    double Avals[5] = {beta, alpha, 1., alpha, beta};

    // long int Nsize = (long int)nx * (long int)ny * (long int)nz;

    // nnz = 5*Nsize;
    // colptr = (int *) malloc ( (Nsize+1) * sizeof(int));
    // rowind = (int *) malloc ( nnz * sizeof(int));
    // nzval = (double *) malloc ( nnz * sizeof(double));

    long int counter = 0;
    colptr[0] = counter;

    for (int iz = 0; iz < nz; ++iz){
        for (int iy = 0; iy < ny; ++iy){
            for (int col = 0; col < nx; ++col){
                for (int j = 0; j < 5; ++j){
                    int row = col + j - 2;
                    rowind[counter] = (row + nx)%nx + iy*nx + iz*nx*ny;
                    nzval [counter] = Avals[j];
                    counter++;
                }
                long int gcol = col + iy*nx + iz*nx*ny;
                colptr[gcol+1] = colptr[gcol]+5;
            }
        }
    }

    // dCreate_CompCol_Matrix (A, Nsize, Nsize, nnz, nzval, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
}

void ddx_rhs(const double * const restrict f, double * restrict df, double dx, int nx, int ny, int nz)
{
    double a10bdx = a10/dx;
    double b10bdx = b10/dx;
    double c10bdx = c10/dx;

    for (int k = 0; k < nz; ++k)
        for (int j = 0; j < ny; ++j)
            for (int i = 0; i < nx; ++i)
                df[LOC(i,j,k)] = ( a10bdx*(f[LOC(i+1,j,k)]-f[LOC(i-1,j,k)]) 
                                 + b10bdx*(f[LOC(i+2,j,k)]-f[LOC(i-2,j,k)]) 
                                 + c10bdx*(f[LOC(i+3,j,k)]-f[LOC(i-3,j,k)]) );
}

double check_error(int nx, int ny, int nz, double *x, double *y, double *z, double *f, double *df)
{
    double err = 0.;
    for (int k = 0; k < nz; ++k){
        for (int j = 0; j < ny; ++j){
            for (int i = 0; i < nx; ++i){
                int idx = LOC(i,j,k);
                err = ( fabs( df[idx] - cos(x[idx]) * cos(y[idx]) * cos(z[idx]) ) > err ) ? fabs( df[idx] - cos(x[idx]) * cos(y[idx]) * cos(z[idx]) ) : err;
            }
        }
    }
    return err;
}

int main(int argc, char *argv[])
{
    
    int nx = 64, ny = 64, nz = 64;

    double dx, dy, dz;
    double *x, *y, *z, *f, *df, *dX;
   
    superlu_vars_t vars;
    
    // int info;

    double *nzval;
    int *rowind, *colptr;
    long int nnz;

    long int Nsize = nx*ny*nz;

    x  = (double *) malloc ( nx * ny * nz * sizeof (double));
    y  = (double *) malloc ( nx * ny * nz * sizeof (double));
    z  = (double *) malloc ( nx * ny * nz * sizeof (double));
    f  = (double *) malloc ( nx * ny * nz * sizeof (double));
    df = (double *) malloc ( nx * ny * nz * sizeof (double));
    dX = (double *) malloc ( nx * ny * nz * sizeof (double));

    dx = 2.*M_PI/nx; dy = 2.*M_PI/ny; dz = 2.*M_PI/nz;

    nnz = 5*Nsize;
    colptr = (int *) malloc ( (Nsize+1) * sizeof(int));
    rowind = (int *) malloc ( nnz * sizeof(int));
    nzval  = (double *) malloc ( nnz * sizeof(double));

    create_pade_matrix(nzval, rowind, colptr, nnz, nx, ny, nz);

    initialize_superlu_vars(nzval, rowind, colptr, Nsize, nnz, df, dX, &vars);

    for (int k = 0; k < nz; ++k){
        for (int j = 0; j < ny; ++j){
            for (int i = 0; i < nx; ++i){
                int idx = LOC(i,j,k);
                x[idx] = i*dx;
                y[idx] = j*dy;
                z[idx] = k*dz;
                f[idx] = sin(x[idx]) * cos(y[idx]) * cos(z[idx]);
            }
        }
    }

    printf("Starting derivative routines\n");

    clock_t start = clock(), diff;
    ddx_rhs(f,dX,dx,nx,ny,nz);
    MatrixSolve(dX, df, nzval, nx, ny, nz, &vars);
    diff = clock() - start;

    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken = %g seconds\n", (double)msec/1000);

    double err = check_error(nx,ny,nz,x,y,z,f,df);

    printf("Maximum error in x derivative = %g\n",err);

    destroy_superlu_vars(&vars);

    free(x );
    free(y );
    free(z );
    free(f );
    free(df);
    free(dX);
    
    // free(nzval);
    // free(rowind);
    // free(colptr);

    return 0;
}
