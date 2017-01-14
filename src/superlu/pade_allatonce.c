#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "slu_ddefs.h"

#define LOC(i,j,k) ( ((i+nx)%nx) + nx*((j+ny)%ny) + nx*ny*((k+nz)%nz) )

#define alpha (1./2. )
#define beta  (1./20.)
#define a10   ((17./12.)  /2.)
#define b10   ((101./150.)/4.)
#define c10   ((1./100.)  /6.)

void create_pade_matrix(double *nzval, int *rowind, int *colptr, long int nnz, int nx, int ny, int nz)
{
    // double *nzval;
    // int *rowind, *colptr;
    // long int nnz;

    double Avals[5] = {beta, alpha, 1., alpha, beta};

    long int Nsize = (long int)nx * (long int)ny * (long int)nz;

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

void ddx_rhs(const double * const f, double *df, double dx, int nx, int ny, int nz)
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

void ddx(const double * const f, double *dX, double *df, double dx, SuperMatrix *A, double *nzval,
         SuperMatrix *L, SuperMatrix *U, SuperMatrix *B, SuperMatrix *X,
         int *perm_c, int *perm_r, int *etree, char *equed, double *R, double *C, void *work, int lwork, double *rpg, double *rcond,
         double *ferr, double *berr, GlobalLU_t *Glu, mem_usage_t *mem_usage, 
         superlu_options_t *options, SuperLUStat_t *stat, int nx, int ny, int nz)
{
    int info;

    ddx_rhs(f,dX,dx,nx,ny,nz);

    ((NCformat*) A->Store)->nzval = (void*)(nzval);

    ((DNformat*) B->Store)->nzval = (void*)(dX);
    ((DNformat*) X->Store)->nzval = (void*)(df);
    
    // dgssv(options, A, perm_c, perm_r, L, U, B, stat, &info);
    dgssvx(options, A, perm_c, perm_r, etree, equed, R, C,
           L, U, work, lwork, B, X, rpg, rcond, ferr, berr,
           Glu, mem_usage, stat, &info);

    // DNformat * tmp = (DNformat *)((*B).Store);
    // df = (double *)((*tmp).nzval);
    // dCreate_Dense_Matrix (B, nx, ny*nz, df, nx, SLU_DN, SLU_D, SLU_GE);
    // dPrint_Dense_Matrix("B", B);
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
    
    char equed[1];
    double *R, *C;
    double *ferr, *berr;

    void *work;

    superlu_options_t options;
    int *perm_c, *perm_r, *etree;
    SuperLUStat_t stat;
    int info, lwork;

    double rpg, rcond;

    GlobalLU_t Glu;
    mem_usage_t mem_usage;

    double *nzval;
    int *rowind, *colptr;
    long int nnz;

    SuperMatrix A, L, U, B, X;

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
    dCreate_CompCol_Matrix (&A, Nsize, Nsize, nnz, nzval, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
    // dPrint_CompCol_Matrix("A", &A);

    perm_r = (int *) malloc( Nsize * sizeof(int));
    perm_c = (int *) malloc( Nsize * sizeof(int));
    etree  = (int *) malloc( Nsize * sizeof(int));

    R = (double *) malloc( A.nrow*sizeof(double) );
    C = (double *) malloc( A.ncol*sizeof(double) );
    ferr = (double *) malloc( 1*sizeof(double) ); // 1 is the no. of RHS
    berr = (double *) malloc( 1*sizeof(double) ); // 1 is the no. of RHS

    lwork = 0;

    set_default_options(&options);
    options.Equil = NO;
    options.DiagPivotThresh = 1.0;
    options.Trans = NOTRANS;
    options.PrintStat = NO;
    options.ColPerm = COLAMD;
    options.RowPerm = NO;

    StatInit(&stat);

    dCreate_Dense_Matrix (&B, Nsize, 1, df, Nsize, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix (&X, Nsize, 1, dX, Nsize, SLU_DN, SLU_D, SLU_GE);
    // dPrint_Dense_Matrix("B", &B);
 
    //dgssv(&options, &A, perm_c, perm_r, &L, &U, &B, &stat, &info);
    
    dgssvx(&options, &A, perm_c, perm_r, etree, equed, R, C,
           &L, &U, work, lwork, &B, &X, &rpg, &rcond, ferr, berr,
           &Glu, &mem_usage, &stat, &info);

    // options.Fact = SamePattern_SameRowPerm;
    options.Fact = FACTORED;
    StatInit(&stat);

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
    ddx(f,dX,df,dx,&A,nzval,
        &L,&U,&B,&X,perm_c,perm_r,etree,equed,R,C,work,lwork,&rpg,&rcond,ferr,berr,&Glu,&mem_usage,&options,&stat,nx,ny,nz);
    diff = clock() - start;

    int msec = diff * 1000 / CLOCKS_PER_SEC;
    printf("Time taken = %g seconds\n", (double)msec/1000);

    double err = check_error(nx,ny,nz,x,y,z,f,df);

    printf("Maximum error in x derivative = %g\n",err);

    Destroy_CompCol_Matrix(&A);
    Destroy_SuperMatrix_Store(&B);
    Destroy_SuperNode_Matrix(&L);
    Destroy_CompCol_Matrix(&U);
    StatFree(&stat);

    free(x );
    free(y );
    free(z );
    free(f );
    free(df);
    free(dX);
    
    free(perm_c);
    free(perm_r);


    return 0;
}
