/*! \file faspcpr.h
 *
 *  \brief Main header file for FASPCPR
 *
 *  \note  This header file contains general constants and data structures of
 *         FASPCPR. It contains macros and data structure definitions; should
 *         not include function declarations here.
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2015--2022 by the FASP team. All rights reserved.
 *---------------------------------------------------------------------------------
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "fasp_const.h"

#ifndef __FASPCPR_HEADER__   /*-- allow multiple inclusions --*/
#define __FASPCPR_HEADER__   /**< indicate faspcpr.h has been included */


/*---------------------------*/
/*---     Parameters      ---*/
/*---------------------------*/

extern INT  last_iter_c;
extern INT  resetup_c;
extern INT  fasp_called_times;
extern REAL total_setup_time; /**< Total setup time                  */
extern REAL total_solve_time; /**< Total setup time                  */


/*---------------------------*/
/*---   Data structures   ---*/
/*---------------------------*/

/**
 * \struct block_STR
 *
 * \brief  Block REAL matrix format for reservoir simulation
 */
typedef struct block_STR {
    
    //! reservoir-reservoir block
    dSTRmat ResRes;
    
    //! reservoir-well block
    dCSRmat ResWel;
    
    //! well-reservoir block
    dCSRmat WelRes;
    
    //! well-well block
    dCSRmat WelWel;
    
} block_STR; /**< Special block matrix for Reservoir Simulation */

/**
 * \struct block_BSR
 *
 * \brief  Block REAL matrix format for reservoir simulation
 */
typedef struct block_BSR {
    
    //! reservoir-reservoir block
    dBSRmat ResRes;
    
    //! reservoir-well block
    dCSRmat ResWel;
    
    //! well-reservoir block
    dCSRmat WelRes;
    
    //! well-well block
    dCSRmat WelWel;
    
} block_BSR; /**< Block of BSR matrices of REAL type */

/**
 * \struct precond_blkoil_data
 *
 * \brief  Data for the preconditioners to reservoir simulation problems in BBSR
 *
 * \note   This is only needed for FIM Black Oil systems with wells
 */
typedef struct {
    
    //! problem data in block_STR format
    block_STR *A;
    
    //! problem data in dBLCmat format
    dBLCmat *Ablc;
    
    //! problem data in CSR format
    dCSRmat *Acsr;
    
    //! level of fill-in for structured ILU(k)
    INT ILU_lfil;
    
    //! LU matrix for Reservoir-Reservoir block in STR format
    dSTRmat *LU;
    
    //! LU matrix for Reservoir-Reservoir block in CSR format
    ILU_data *LUcsr;
    
    //! AMG data for presure-presure block
    AMG_data *mgl_data;
    
    //! print level in AMG preconditioner
    SHORT print_level;
    
    //! max number of iterations of AMG preconditioner
    INT maxit_AMG;
    
    //! max number of AMG levels
    SHORT max_levels;
    
    //! tolerance for AMG preconditioner
    REAL amg_tol;
    
    //! AMG cycle type
    SHORT cycle_type;
    
    //! AMG smoother type
    SHORT smoother;
    
    //! number of presmoothing
    SHORT presmooth_iter;
    
    //! number of postsmoothing
    SHORT postsmooth_iter;
    
    //! coarsening type
    SHORT coarsening_type;
    
    //! relaxation parameter for SOR smoother
    REAL relaxation;
    
    //! switch of scaling of coarse grid correction
    SHORT coarse_scaling;
    
    //! max number of iterations
    INT  maxit;
    
    //! number of iterations for restart
    INT restart;
    
    //! tolerance for convergence
    REAL tol;
    
    //! inverse of the Schur complement (-I - Awr*Arr^{-1}*Arw)^{-1}
    REAL *invS; // in which Arr may be replaced by LU
    
    //! Diag(PS) * inv(Diag(SS))
    dvector *DPSinvDSS;
    
    //! whether the matirx is scaled or not
    SHORT    scaled;
    
    //! the diagonal for diagonal scaling */
    precond_diag_str *diag;
    
    //! the inverse of the diagonal for GS/block GS smoother (whole reservoir)
    dvector *diaginv;
    
    //! the pivot for the GS/block GS smoother (whole reservoir)
    ivector *pivot;
    
    //! the inverse of the diagonals for GS/block GS smoother (saturation)
    dvector *diaginvS;
    
    //! the pivot for the GS/block GS smoother (saturation)
    ivector *pivotS;
    
    //! order for smoothing
    ivector *order;
    
    //! index for perforations
    ivector *perf_idx;
    
    // block structure of coefficient matrix
    dSTRmat *RR; /**< Diagonal scaled reservoir block */
    dCSRmat *WW; /**< Argumented well block */
    dCSRmat *PP; /**< pressure block after diagonal scaling */
    dSTRmat *SS; /**< saturation block after diaogonal scaling */
    
    // temporary work space
    dvector r;   /**< temporary dvector used to store and restore the residual */
    REAL   *w;   /**< temporary work space for other usage */
    
} precond_blkoil_data; /**< Precond data for Reservoir Simulation */

/**
 * \struct precond_blkoil_bbsr_data
 *
 * \brief  Data for the preconditioners to reservoir simulation problems in BBSR
 *
 * \note   This is only needed for the Black Oil model with wells
 */
typedef struct {
    
    //-------------------------------------------------------------------------------
    //! Part 1: Basic data
    //-------------------------------------------------------------------------------
    block_BSR *A; /**< whole jacobian system in block_BSRmat */
    
    //-------------------------------------------------------------------------------
    //! Part 2: Data for CPR-like preconditioner for reservoir block
    //-------------------------------------------------------------------------------
    // diagonal scaling for reservoir block
    SHORT scaled; /**< scaled = 1 means the RR block is diagonal scaled */
    dBSRmat *RR;  /**< reservoir block */
    dvector *diaginv_noscale; /**< inverse of block diagonal for scaling */
    
    // neighborhood and reordering of reservoir block
    ivector *neigh;     /**< neighbor information of the reservoir block */
    ivector *order;     /**< ordering of the reservoir block */
    
    // data for GS/bGS smoother for saturation block
    dBSRmat *SS;        /**< saturation block */
    dvector *diaginv_S; /**< inverse of the diagonal blocks of saturation block */
    ivector *pivot_S;   /**< pivoting for the GS smoothers for saturation block */
    
    // data of ILU for saturation block
    ILU_data *LU_S;     /**< ILU setup data for saturation block */
    
    // data of AMG for pressure block
    dCSRmat *PP;        /**< pressure block */
    AMG_data *mgl_data; /**< AMG data for presure-presure block */
    
    //! data of ILU setup for pressure block
    ILU_data *LU_P;
    
    //! print level in AMG preconditioner
    SHORT print_level;
    
    //! max number of iterations of AMG preconditioner
    INT maxit_AMG;
    
    //! max number of AMG levels
    SHORT max_levels;
    
    //! tolerance for AMG preconditioner
    REAL amg_tol;
    
    //! AMG cycle type
    SHORT cycle_type;
    
    //! AMG smoother type
    SHORT smoother;
    
    //! AMG smoothing order
    SHORT smooth_order;
    
    //! number of presmoothing
    SHORT presmooth_iter;
    
    //! number of postsmoothing
    SHORT postsmooth_iter;
    
    //! coarsening type
    SHORT coarsening_type;
    
    //! coarset dof
    INT coarse_dof;
    
    //! coarse level solver type
    SHORT coarse_solver;
    
    //! relaxation parameter for SOR smoother
    REAL relaxation;
    
    //! switch of scaling of coarse grid correction
    SHORT coarse_scaling;
    
    //! degree of the polynomial used by AMLI cycle
    SHORT amli_degree;
    
    //! coefficients of the polynomial used by AMLI cycle
    REAL *amli_coef;
    
    //! relaxation parameter for smoothing the tentative prolongation
    REAL tentative_smooth;
    
    // data of GS/bGS smoother for reservoir block
    dvector *diaginv; /**< inverse of the diagonal blocks of reservoir block */
    ivector *pivot;   /**< pivot for the GS smoothers for the reservoir matrix */
    
    //! data of ILU for reservoir block
    ILU_data *LU;
    
    // data for the argumented well block
    //! index of blocks which have perforation
    ivector *perf_idx;
    
    //! index of blocks which are neighbors of perforations (include perforations)
    ivector *perf_neigh;
    
    //! Argumented well block
    dCSRmat *WW;
    
    //! data for direct solver for argumented well block
    void *Numeric;
    
    //! inverse of the schur complement (-I - Awr*Arr^{-1}*Arw)^{-1}
    REAL *invS;  // in which Arr may be replaced by LU
    
    // parameters for krylov method used for blocks
    INT  maxit;    /**< max number of iterations */
    INT  restart;  /**< number of iterations for restart */
    REAL tol;      /**< tolerance */
    
    // temporary work space
    dvector r; /**< temporary dvector used to store and restore the residual */
    REAL *w;   /**< temporary work space for other usage */
    
} precond_blkoil_bbsr_data; /**< FASP precond data for Reservoir Simulation */

#endif /* end if for __FASPCPR_HEADER__ */

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
