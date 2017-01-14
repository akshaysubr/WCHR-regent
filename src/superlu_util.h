#ifndef __SUPERLU_UTIL__
#define __SUPERLU_UTIL__

#include "slu_ddefs.h"

typedef struct superlu_vars_t {
    SuperMatrix A;
    SuperMatrix L;
    SuperMatrix U;
    SuperMatrix B;
    SuperMatrix X;
    GlobalLU_t Glu;
    mem_usage_t mem_usage;
    superlu_options_t options;
    SuperLUStat_t stat;

    int *perm_c;
    int *perm_r;
    int *etree;

    double *R;
    double *C;
    double *ferr;
    double *berr;

    double rpg;
    double rcond;

    void *work;
    int  lwork;

    char equed[1];
} superlu_vars_t;

void MatrixSolve(double * restrict dX, double * restrict df, double * restrict nzval, int nx, int ny, int nz, superlu_vars_t *vars);

void destroy_superlu_vars(superlu_vars_t *vars);

superlu_vars_t initialize_superlu_vars(double *nzval, int* rowind, int *colptr, long int Nsize, long int nnz, double *b, double *x);

#endif
