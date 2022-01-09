/*! \file  PreBSR.c
 *
 *  \brief Preconditioners for petroleum reservoir simulation in BSR format
 *
 *  \note  The routines in this file are designed for the reservoir block. Implicit
 *         wells should not exist or they have been combined with the pressure
 *         unknowns by introducing auxiliary saturation variables!
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2015--2022 by the FASP team. All rights reserved.
 *---------------------------------------------------------------------------------
 */

#ifdef _OPENMP
#include <omp.h>
#endif

#include "fasp.h"
#include "fasp_functs.h"
#include "faspcpr.h"
#include "faspcpr_functs.h"

/*---------------------------------*/
/*--   Decoration of Functions   --*/
/*---------------------------------*/

static void compute_rp (REAL *, dBSRmat *, REAL *, REAL *);
void fasp_aux_pres2blkoil (INT , REAL *, INT , REAL *, INT , INT );


/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn void fasp_precond_dbsr_CPR (REAL *r, REAL *z, void *data)
 *
 * \brief Get z = Br by FASP preconditioner for dBSRmat matrices -- general
 *
 * \param r      Pointer to residual
 * \param z      Pointer to preconditioned residual
 * \param data   Pointer to precondition data
 *
 * \author Li Zhao
 * \date   03/19/2021
 */
void fasp_precond_dbsr_CPR (REAL *r,
                            REAL *z,
                            void *data)
{
    precond_blkoil_bbsr_data *predata = (precond_blkoil_bbsr_data *)data;
    dBSRmat  *A     = predata->RR;
    AMG_data *mgl   = predata->mgl_data;
    ivector  *order = predata->order;
    ILU_data *LU    = predata->LU;
    
    const INT maxit = predata->maxit; // it refers to AMG_maxit here, Zhiyang Zhou
    const INT nc    = A->nb;
    const INT ngrid = A->ROW;
    const INT n     = nc*ngrid;     // whole size
    const INT m     = mgl[0].A.row; // pressure block size
    
    INT i;
    
    dvector rr; rr.row = n; rr.val = r;
    dvector zz; zz.row = n; zz.val = z;
    
    memset(z, 0X0, n*sizeof(REAL));
    
    //-----------------------------------------------------
    // Stage 2: solve pressure block
    //-----------------------------------------------------
    
    // 2.1: Restrict tempr to rp
    compute_rp(r, A, z, mgl->b.val);
    
    // 2.2: solve pressure block
    AMG_param amgparam;
    fasp_param_amg_init(&amgparam);
    amgparam.cycle_type      = predata->cycle_type;
    amgparam.coarse_solver   = predata->coarse_solver;
    amgparam.smoother        = predata->smoother;
    amgparam.presmooth_iter  = predata->presmooth_iter;
    amgparam.postsmooth_iter = predata->postsmooth_iter;
    amgparam.relaxation      = predata->relaxation;
    amgparam.coarse_scaling  = predata->coarse_scaling;
    amgparam.ILU_levels      = predata->mgl_data->ILU_levels;
    
    memset(mgl->x.val, 0X0, m*sizeof(REAL));
    
    switch (amgparam.cycle_type)
    {
        case AMLI_CYCLE:
            for (i=0;i<maxit;++i) fasp_solver_amli(mgl,&amgparam,0);
            break;
        case NL_AMLI_CYCLE:
            for (i=0; i<maxit; ++i) fasp_solver_namli(mgl,&amgparam,0,mgl[0].num_levels);
            break;
        default:
            for (i = 0; i <maxit; ++i) fasp_solver_mgcycle(mgl, &amgparam);
            break;
    }
    
    // 2.3: Prolongate zp to z
    fasp_aux_pres2blkoil(m, mgl->x.val, n, z, nc, 1);
    
    //-----------------------------------------------------
    // Stage 3: ILU for whole matrix
    //-----------------------------------------------------
    fasp_smoother_dbsr_ilu(A, &rr, &zz, LU);
}



/**
 * \fn void fasp_aux_pres2blkoil (INT Npres, REAL *pres, INT Nx, REAL *x,
 *                                INT jneq, INT order)
 *
 * \brief Prolongate pressure vector to whole vector
 *
 * \param Npres  length of the pressure vector
 * \param pres   pressure vector
 * \param Nx     length of whole vector
 * \param x      whole vector: pressure should always be the first variables
 * \param jneq   number of components
 * \param order  ordering used for the whole matrix
 *
 * \author Shiquan Zhang
 * \date 07/16/2010
 *
 * \author Feng Chunsheng, Yue Xiaoqiang on 03/10/2011
 */
void fasp_aux_pres2blkoil (INT Npres,
                           REAL *pres,
                           INT Nx,
                           REAL *x,
                           INT jneq,
                           INT order)
{
    INT i;
    
    switch ( order ) {
            
        case 1: // grid-wise ordering
            
#ifdef _OPENMP
            if (Npres > OPENMP_HOLDS) {
#pragma omp parallel for private(i)
                for (i=0; i<Npres; ++i) x[i*jneq]=pres[i];
            }
            else {
#endif
                for (i=0; i<Npres; ++i) x[i*jneq]=pres[i];
#ifdef _OPENMP
            }
#endif
            break;
            
        case 2: // component-wise ordering
            
#ifdef _OPENMP
            if (Npres > OPENMP_HOLDS) {
#pragma omp parallel for private(i)
                for (i=0; i<Npres; ++i) x[i]=pres[i];
            }
            else {
#endif
                for (i=0; i<Npres; ++i) x[i]=pres[i];
#ifdef _OPENMP
            }
#endif
            break;
            
        default: // unknown ordering
            
            printf("### ERROR: Unknown ordering type %d for Blkoil!\n", order);
            exit(ERROR_INPUT_PAR);
            break;
    }
}


/*---------------------------------*/
/*--      Private Functions      --*/
/*---------------------------------*/

/**
 * \fn static void get_temp (INT  nc, REAL *val, REAL *z, INT pos, INT zrow)
 *
 * \brief Compute temp used in the subroutine compute_rp
 *
 * \param nc    number of components
 * \parma val   pointer to val
 * \param z     pointer to z
 * \param pos   true position in the value array (i-th block * block size)
 * \param zrow  true column index in the z array (i-th column * number of components)
 *
 * \author Xiaozhe Hu
 * \date   12/25/2010
 *
 * \note work for general nb (Xiaozhe)
 * \note modified by Chensong on 11/09/2011
 */
static REAL get_temp (INT  nc,
                      REAL *val,
                      REAL *z,
                      INT pos,
                      INT zrow)
{
    REAL temp = 0.0;
    
    switch (nc){
            
        case 2:
            temp += val[pos+1]*z[zrow+1];
            break;
            
        case 3:
            temp += val[pos+1]*z[zrow+1] + val[pos+2]*z[zrow+2];
            break;
            
        case 5:
            temp += val[pos+1]*z[zrow+1] + val[pos+2]*z[zrow+2] + \
            val[pos+3]*z[zrow+3] + val[pos+4]*z[zrow+4];
            break;
            
        case 7:
            temp += val[pos+1]*z[zrow+1] + val[pos+2]*z[zrow+2] + \
            val[pos+3]*z[zrow+3] + val[pos+4]*z[zrow+4] + \
            val[pos+5]*z[zrow+5] + val[pos+6]*z[zrow+6];
            break;
            
        default:
        {
            INT k;
            for (k = 1; k < nc; k ++) temp += val[pos+k]*z[zrow+k];
        }
            break;
            
    }
    
    return temp;
}


/**
 * \fn static void compute_rp (REAL *r, dBSRmat *A, REAL *z, REAL *rp)
 *
 * \brief Compute rp := r_p - A_{ps}*z_s, where r_p is the pressure part of r,
 *        and A_{ps} is the top right part of A.
 *
 * \param r  pointer to the whole residual
 * \parma A  pointer to the whole coefficient matrix
 * \param z  pointer to the whole current solution
 * \param rp pointer to resulting REAL array
 *
 * \author Xiaozhe Hu, Zhiyang Zhou
 * \date   12/25/2010
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/29/2012
 */
static void compute_rp (REAL *r,
                        dBSRmat *A,
                        REAL *z,
                        REAL *rp)
{
    const INT nc       = A->nb;
    const INT ngrid    = A->ROW;
    const INT jump     = nc*nc;
    
    INT    *IA       = A->IA;
    INT    *JA       = A->JA;
    REAL *val      = A->val;
    REAL  temp     = 0.0;
    
    INT i, j, ibegin, iend;
    INT col, zrow, pos;
    
#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend;
    INT nthreads = fasp_get_num_threads();
#endif
    
#ifdef _OPENMP
    if (ngrid > OPENMP_HOLDS) {
#pragma omp parallel for private(myid,mybegin,myend,i,j,ibegin,iend,col,pos,zrow,temp)
        for (myid = 0; myid < nthreads; myid++) {
            fasp_get_start_end(myid, nthreads, ngrid, &mybegin, &myend);
            for (i = mybegin; i < myend; i++) {
                temp = 0.0;
                ibegin = IA[i]; iend = IA[i+1];
                for (j = ibegin; j < iend; j ++) {
                    col  = JA[j];
                    pos  = j*jump;
                    zrow = col*nc;
                    temp += get_temp(nc, val, z, pos, zrow);
                }
                rp[i] = r[i*nc] - temp;
            }
        }
    }
    else {
#endif
        for (i = 0; i < ngrid; ++i) {
            temp = 0.0;
            ibegin = IA[i]; iend = IA[i+1];
            for (j = ibegin; j < iend; j ++) {
                col  = JA[j];
                pos  = j*jump;
                zrow = col*nc;
                temp += get_temp(nc, val, z, pos, zrow);
            }
            rp[i] = r[i*nc] - temp;
        }
#ifdef _OPENMP
    }
#endif
}


/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
