#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#include "slu_ddefs.h"
#include "superlu_util.h"

#define LOC(i,j,k) ( ((i+nx)%nx) + nx*((j+ny)%ny) + nx*ny*((k+nz)%nz) )

void MatrixSolve(double * restrict dX, double * restrict df, double * restrict nzval, int nx, int ny, int nz, superlu_vars_t *vars)
{
    int info;

    ((NRformat*) vars->A.Store)->nzval = (void*)(nzval);

    ((DNformat*) vars->B.Store)->nzval = (void*)(dX);
    ((DNformat*) vars->X.Store)->nzval = (void*)(df);
    
    // dPrint_CompCol_Matrix("A: ", &(vars->A));
    dgssvx(&(vars->options), &(vars->A), vars->perm_c, vars->perm_r, vars->etree, vars->equed, vars->R, vars->C,
           &(vars->L), &(vars->U), vars->work, vars->lwork, &(vars->B), &(vars->X), &(vars->rpg), &(vars->rcond), vars->ferr, vars->berr,
           &(vars->Glu), &(vars->mem_usage), &(vars->stat), &info);

}

void destroy_superlu_vars(superlu_vars_t *vars)
{
    // Destroy_CompCol_Matrix(&(vars->A));
    Destroy_CompRow_Matrix(&(vars->A));
    Destroy_SuperMatrix_Store(&(vars->B));
    Destroy_SuperMatrix_Store(&(vars->X));
    Destroy_SuperNode_Matrix(&(vars->L));
    Destroy_CompCol_Matrix(&(vars->U));
    StatFree(&(vars->stat));

    SUPERLU_FREE(vars->perm_c);
    SUPERLU_FREE(vars->perm_r);
    SUPERLU_FREE(vars->etree);

    SUPERLU_FREE(vars->R);
    SUPERLU_FREE(vars->C);
    SUPERLU_FREE(vars->ferr);
    SUPERLU_FREE(vars->berr);

    if ( (vars->lwork) != 0 )
        SUPERLU_FREE(vars->work);
}


void initialize_superlu_vars(double *nzval, int* rowind, int *colptr, long int Nsize, long int nnz, double *b, double *x, superlu_vars_t *vars)
{
    int info;

    // dCreate_CompCol_Matrix (&(vars->A), Nsize, Nsize, nnz, nzval, rowind, colptr, SLU_NC, SLU_D, SLU_GE);
    dCreate_CompRow_Matrix (&(vars->A), Nsize, Nsize, nnz, nzval, rowind, colptr, SLU_NR, SLU_D, SLU_GE);
    // dPrint_CompCol_Matrix("A: ", &(vars->A));

    vars->perm_r = (int *) malloc( Nsize * sizeof(int));
    vars->perm_c = (int *) malloc( Nsize * sizeof(int));
    vars->etree  = (int *) malloc( Nsize * sizeof(int));

    vars->R    = (double *) malloc( vars->A.nrow*sizeof(double) );
    vars->C    = (double *) malloc( vars->A.ncol*sizeof(double) );
    vars->ferr = (double *) malloc( 1*sizeof(double) ); // 1 is the no. of RHS
    vars->berr = (double *) malloc( 1*sizeof(double) ); // 1 is the no. of RHS

    vars->lwork = 0;

    set_default_options(&(vars->options));
    vars->options.Equil = YES;
    vars->options.DiagPivotThresh = 1.0;
    vars->options.Trans = NOTRANS;
    vars->options.PrintStat = NO;
    vars->options.ColPerm = COLAMD;
    vars->options.RowPerm = NO;
    vars->options.IterRefine = NOREFINE;
    // vars->options.IterRefine = SLU_EXTRA;

    StatInit(&(vars->stat));

    dCreate_Dense_Matrix (&(vars->B), Nsize, 1, b, Nsize, SLU_DN, SLU_D, SLU_GE);
    dCreate_Dense_Matrix (&(vars->X), Nsize, 1, x, Nsize, SLU_DN, SLU_D, SLU_GE);
 
    dgssvx(&(vars->options), &(vars->A), vars->perm_c, vars->perm_r, vars->etree, vars->equed, vars->R, vars->C,
           &(vars->L), &(vars->U), vars->work, vars->lwork, &(vars->B), &(vars->X), &(vars->rpg), &(vars->rcond), vars->ferr, vars->berr,
           &(vars->Glu), &(vars->mem_usage), &(vars->stat), &info);

    vars->options.Fact = SamePattern_SameRowPerm;
    // vars->options.Fact = FACTORED;

}
