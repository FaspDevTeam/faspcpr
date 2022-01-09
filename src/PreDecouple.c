/*! \file  PreDecouple.c
 *
 *  \brief General decoupling methods for Black Oil FIM Jacobian matrices
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2015--2022 by the FASP team. All rights reserved.
 *---------------------------------------------------------------------------------
 */

#include <math.h>

#include "fasp.h"
#include "fasp_functs.h"

/*---------------------------------*/
/*--  Declare Private Functions  --*/
/*---------------------------------*/

static inline void fasp_blas_sarry_xy_nc2 (const REAL *,const REAL *,REAL *);
static inline void fasp_blas_sarry_xy_nc3 (const REAL *,const REAL *,REAL *);
static inline void smat_solution_nc2      (const REAL *,REAL *);
static inline void smat_solution_nc3      (const REAL *,REAL *);
static inline void fasp_blas_dbsr_col_accumulation (dBSRmat *, REAL *);


/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dBSRmat fasp_dbsr_timpes_decoupling2 (dBSRmat *A, REAL *diaginv, REAL *sol)
 *
 * \brief Compute True-IMPES decoupling matrix.
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param diaginv  Pointer to the True-IMPES decoupling matrix
 * \param sol      Result after decoupling step
 *
 * \return BSR matrix after true impes decoupling
 *
 * \author Zheng Li, Chensong Zhang
 * \date   06/12/2014
 */
dBSRmat fasp_dbsr_timpes_decoupling2 (dBSRmat *A,
                                      REAL *diaginv,
                                      REAL *sol)
{
    dBSRmat B;
    // members of A
    INT     ROW = A->ROW;
    INT     ROW_plus_one = ROW+1;
    INT     COL = A->COL;
    INT     NNZ = A->NNZ;
    INT     nb  = A->nb;
    INT    *IA  = A->IA;
    INT    *JA  = A->JA;
    REAL   *val = A->val;
    
    INT    *IAb  = NULL;
    INT    *JAb  = NULL;
    REAL   *valb = NULL;
    
    INT nb2  = nb*nb;
    INT i,j,k,m;
    
    // Create a dBSRmat 'B'
    B = fasp_dbsr_create(ROW, COL, NNZ, nb, 0);
    
    IAb  = B.IA;
    JAb  = B.JA;
    valb = B.val;
    
    fasp_iarray_cp(ROW_plus_one, IA, IAb);
    fasp_iarray_cp(NNZ, JA, JAb);
    fasp_darray_set(NNZ*nb2, valb, 0.0);
    
    // get A^
    fasp_darray_set(ROW*nb2, diaginv, 0.0);
    fasp_blas_dbsr_col_accumulation(A, diaginv);
    
    switch (nb){
            
        case 2:
            for (i = 0; i < ROW; ++i) {
                smat_solution_nc2(diaginv+i*4, sol+i);
                for (k = IA[i]; k < IA[i+1]; ++k) {
                    fasp_blas_sarry_xy_nc2(val+k*4, sol+i, valb+k*4);
                }
            }// end of main loop
            
            break;
            
        case 3:
            for (i = 0; i < ROW; ++i) {
                smat_solution_nc3(diaginv+i*9, sol+i*2);
                for (k = IA[i]; k < IA[i+1]; ++k) {
                    fasp_blas_sarry_xy_nc3(val+k*9, sol+i*2, valb+k*9);
                }
            }// end of main loop
            
            break;
    }
    
    return B;
}


/**
 * \fn void fasp_precond_dbsr_timpes_decoupling2 (REAL *val, REAL *sol,
 *                                                INT n, INT nb, REAL *pb)
 *
 * \brief Compute first row of z:= Br, where B is True-IMPES decoupling matrix.
 *
 * \param val    Pointer to origin right hand side.
 * \param sol    Pointer to array formed by True-IMPES decoupling.
 * \param n      Number of block rows.
 * \param nb     Number of components.
 * \param pb     Pointer to decoupled right hand side.
 *
 * \author Zheng Li, Chensong Zhang
 * \date   06/12/2014
 *
 * \note   Work for general nb
 */
void fasp_precond_dbsr_timpes_decoupling2 (REAL *val,
                                           REAL *sol,
                                           INT n,
                                           INT nb,
                                           REAL *pb)
{
    INT i;
    
    switch (nb) {
        case 2:
            for (i=0; i<n; i++) {
                pb[2*i]   = val[2*i] + sol[i]*val[2*i+1];
                pb[2*i+1] = val[2*i+1];
            }
            break;
            
        case 3:
            for (i=0; i<n; i++) {
                pb[3*i]   = val[3*i] + sol[i*2]*val[3*i+1] + sol[i*2+1]*val[3*i+2];
                pb[3*i+1] = val[3*i+1];
                pb[3*i+2] = val[3*i+2];
            }
            
            break;
    }
}


/*---------------------------------*/
/*--     Private Functions       --*/
/*---------------------------------*/

/**
 * \fn static inline void fasp_blas_sarry_xy_nc2 (const REAL *val, const REAL *sol, REAL *bval)
 *
 * \brief bval = sol*val.
 *
 * \param val     Pointer to the small matrix
 * \param sol     Pointer to solution
 *
 * \author Zheng Li, Chensong Zhang
 * \date   06/12/2014
 */
static inline void fasp_blas_sarry_xy_nc2 (const REAL *val,
                                           const REAL *sol,
                                           REAL *bval)
{
    const REAL a11 = val[0], a12 = val[1];
    const REAL a21 = val[2], a22 = val[3];
    
    bval[0] = a11 + sol[0]*a21;
    bval[1] = a12 + sol[0]*a22;
    bval[2] = a21;
    bval[3] = a22;
}


/**
 * \fn static inline void smat_solution_nc2 (const REAL *val, const REAL *sol, REAL *bval)
 *
 * \brief bval = sol*val.
 *
 * \param val     Pointer to the small matrix
 * \param sol     Pointer to solution
 *
 * \author Zheng Li, Chensong Zhang
 * \date   06/12/2014
 */
static inline void fasp_blas_sarry_xy_nc3 (const REAL *val,
                                           const REAL *sol,
                                           REAL *bval)
{
    const REAL a11 = val[0], a12 = val[1], a13 = val[2];
    const REAL a21 = val[3], a22 = val[4], a23 = val[5];
    const REAL a31 = val[6], a32 = val[7], a33 = val[8];
    
    const REAL s0 = sol[0], s1 = sol[1];
    
    bval[0] = a11 + s0*a21 + s1*a31;
    bval[1] = a12 + s0*a22 + s1*a32;
    bval[2] = a13 + s0*a23 + s1*a33;
    bval[3] = a21;
    bval[4] = a22;
    bval[5] = a23;
    bval[6] = a31;
    bval[7] = a32;
    bval[8] = a33;
}


/**
 * \fn static inline void smat_solution_nc2 (const REAL *val, REAL *sol)
 *
 * \brief a12 + sol[0]*a22 = 0.
 *
 * \param val     Pointer to the small matrix
 * \param sol     Pointer to solution
 *
 * \author Zheng Li, Chensong Zhang
 * \date   06/12/2014
 */
static inline void smat_solution_nc2 (const REAL *val,
                                      REAL *sol)
{
    REAL a11 = val[0], a12 = val[1];
    REAL a21 = val[2], a22 = val[3];
    
    if (ABS(a22) < SMALLREAL) a22 = SMALLREAL;
    sol[0] = -a12/a22;
}

/**
 * \fn static inline void smat_solution_nc3 (const REAL *val, REAL *sol)
 *
 * \brief a12+x1*a22+x2*a32 = 0, a13+x1*a23+x2*a33 = 0.
 *
 * \param val     Pointer to the small matrix
 * \param sol     Pointer to solution
 *
 * \author Zheng Li, Chensong Zhang
 * \date   06/12/2014
 */
static inline void smat_solution_nc3 (const REAL *val,
                                      REAL *sol)
{
    const REAL a11 = val[0], a12 = val[1], a13 = val[2];
    const REAL a21 = val[3], a22 = val[4], a23 = val[5];
    const REAL a31 = val[6], a32 = val[7], a33 = val[8];
    
    REAL det  = a22*a33 - a23*a32;
    REAL tmp1 = a13*a32 - a12*a33;
    REAL tmp2 = a12*a23 - a22*a13;
    
    if (ABS(det) < SMALLREAL)  det = SMALLREAL;
    
    sol[0] = tmp1/det;
    sol[1] = tmp2/det;
}


/**
 * \fn static inline void fasp_blas_dbsr_col_accumulation (dBSRmat *A, REAL *bval)
 *
 * \brief Compute the sum of cloumn blocks.
 *
 * \param A        Pointer to the dBSRmat matrix
 * \param bval     Pointer to the sum of cloumn blocks
 *
 * \author Zheng Li, Chensong Zhang
 * \date   06/12/2014
 */
static inline void fasp_blas_dbsr_col_accumulation (dBSRmat *A,
                                                    REAL *bval)
{
    INT *ia = A->IA;
    INT *ja = A->JA;
    REAL *val = A->val;
    INT row = A->ROW;
    INT nb = A->nb;
    INT i,j,jj,kk;
    
    switch (nb) {
        case 2:
            for(i=0; i<row; ++i){
                for (j=ia[i]; j<ia[i+1];++j) {
                    kk = ja[j]*4;
                    jj = j*4;
                    bval[kk]   += val[jj],   bval[kk+1] += val[jj+1];
                    bval[kk+2] += val[jj+2], bval[kk+3] += val[jj+3];
                }
            }
            break;
            
        case 3:
            for(i=0; i<row; ++i){
                for (j=ia[i]; j<ia[i+1];++j) {
                    kk = ja[j]*9;
                    jj = j*9;
                    bval[kk]   += val[jj],   bval[kk+1] += val[jj+1], bval[kk+2] += val[jj+2];
                    bval[kk+3] += val[jj+3], bval[kk+4] += val[jj+4], bval[kk+5] += val[jj+5];
                    bval[kk+6] += val[jj+6], bval[kk+7] += val[jj+7], bval[kk+8] += val[jj+8];
                }
            }
            break;
    }
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
