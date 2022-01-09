/*! \file  fasp4blkoil.c
 *
 *  \brief C/C++ Interface of FASP solvers for BlkOil simulator
 *
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2015--2020 by the FASP team. All rights reserved.
 *---------------------------------------------------------------------------------
 *
 *  \warning Do NOT use auto-indentation for this file!
 */

#include <time.h>


#include "fasp.h"
#include "fasp_block.h"
#include "fasp_functs.h"
#include "faspcpr.h"
#include "faspcpr_functs.h"


/*---------------------------------*/
/*--      Default Setting        --*/
/*---------------------------------*/
REAL total_linear_time = 0.0; /**< Total used time of linear solvers */
REAL total_start_time  = 0.0; /**< Total used time                   */
REAL total_setup_time  = 0.0; /**< Total setup time                  */
REAL total_solve_time  = 0.0; /**< Total setup time                  */
INT  total_iter        = 0;   /**< Total number of iterations        */
INT  fasp_called_times = 0;   /**< Tatal FASP calls                  */

INT  last_iter_c   = 0;
INT  resetup_c     = 1;
INT  last_row_c    = 0;
INT  last_nnz_c    = 0;


/**
 * \fn void FASP_BSRSOL (INT *nrowb, INT *ncolb, INT *nnzb, INT *ia, INT *ja,
 *                       REAL *a, INT *nb, INT *nrow, REAL *b, REAL *x,
 *                       INT *maxit, REAL *tol, INT *itmethod, INT *restart,
 *                       INT *fillin, INT *output, INT *noconv)
 *
 * \brief The C interface to FASP solvers in dBSRmat matrix format
 *
 * \param nrowb       Number of row blocks of coefficient matrix
 * \param ncolb       Number of column blocks of coefficient matrix
 * \param nnzb        Number of nonzero blocks of coefficient matrix
 * \param ia          Row indices for dBSRmat matrix A
 * \param ja          Column positions for dBSRmat matrix A
 * \param a           Value of nonzeros of dBSRmat matrix A
 * \param nb          Number of blocks
 * \param nrow        Number of rows of A
 * \param b           Right-hand side vector
 * \param x           Solution vector (IN: Initial guess, OUT: approx. solution)
 * \param maxit       Max number of iterations
 * \param tol         Convergence tolerance
 * \param itmethod    Iterative method type
 * \param restart     Restarting number for GMRES-type methods
 * \param fillin      Level of fill-in's for ILU when ILU is used
 * \param output      Print level for output
 * \param noconv      Whether the method converges
 *
 * \author Chensong Zhang, Li Zhao
 * \date   11/09/2021
 * 
 */
void FASP_BSRSOL_ASCPR (dBSRmat   *A,
                        dvector   *b,
                        dvector   *x,
                        ITS_param *itparam,
                        ILU_param *iluparam,
                        AMG_param *amgparam,
                        const INT mu)
{ 
    INT iter;
    INT print_level = itparam->print_level;
    
   if ( last_row_c != A->ROW ) {
        resetup_c = 1;
        printf("last row %d %d, last nnz %d %d\n", last_row_c, A->ROW,last_nnz_c , A->NNZ );
        last_row_c = A->ROW;
        last_nnz_c = A->NNZ;
    }

    if ( print_level > 9 ) {
        fasp_dbsr_write("A_bsr.dat",A);
        fasp_dvec_write("b.dat",b);
    }
    

    /*
    *  True-IMPES decoupling method
    */
    // diagonal scaling
    dvector diaginv;
    fasp_dvec_alloc(A->ROW*A->nb*A->nb, &diaginv);
    
    // True-IMPES decoupling
    dvector fsc;
    fasp_dvec_alloc(b->row, &fsc);

    REAL *coeff = (REAL*) fasp_mem_calloc(A->ROW*(A->nb-1), sizeof(REAL));
    dBSRmat Asc = fasp_dbsr_timpes_decoupling2(A, diaginv.val, coeff);
    fasp_precond_dbsr_timpes_decoupling2(b->val, coeff, A->ROW, A->nb, fsc.val);
    fasp_dvec_free(&diaginv);
    fasp_mem_free(coeff);


    ivector order;
    fasp_ivec_alloc(A->ROW, &order);
    
    // ordering
    INT i;
        
#ifdef _OPENMP
#pragma omp parallel for private(i)
#endif
    for ( i = 0; i < A->ROW; ++i ) order.val[i] = i;
        

	
    if (!fasp_called_times) printf("threshold mu = %d\n", mu); 
    iter = fasp_solver_dbsr_krylov_ASCPR (&Asc, &fsc, x, itparam, iluparam, amgparam, NULL, &order);
                                    

    // Clean up   
    fasp_dbsr_free(&Asc);
    fasp_dvec_free(&fsc);
    fasp_ivec_free(&order);

    // Record the number of last iterations
    if ( iter > 0 ) last_iter_c= iter;
    else last_iter_c = itparam->maxit; 


    if (last_iter_c > mu ) resetup_c = 1;
    else{
 	    resetup_c = 0;
    }

    if ( print_level > 9 ) {
		fasp_dvec_write("sol.dat", x);
        printf("Press ENTER to continue..."); getchar();
    }

FINISHED:

    ++fasp_called_times;       
    
    return;
}

 
/**
 * \fn INT fasp_solver_dbsr_krylov_ASCPR (dBSRmat *A, dvector *b, dvector *x,
*                                         ITS_param *itparam, ILU_param *iluparam,
*                                         AMG_param *amgparam, ivector *neigh,
*                                         ivector *order)
*
 * \brief Solve Ax=b by ASCPR preconditioned Krylov methods for BSR matrices
 *        in fully-implicit reservoir simulation
 *
 * \param A         Pointer to the coeff matrix in dBSRmat format
 * \param b         Pointer to the right hand side in dvector format
 * \param x         Pointer to the approx solution in dvector format
 * \param itparam   Pointer to parameters for iterative solvers
 * \param iluparam  Pointer to parameters for ILU
 * \param amgparam  pointer to parameters for AMG
 * \param neigh     Pointer to the information of neighborhood
 * \param order     Pointer to the specified smoothing order
 *
 * \return          Number of iterations if succeed
 *
 * \author Chensong Zhang, Li Zhao
 * \date   11/09/2021
 */
INT fasp_solver_dbsr_krylov_ASCPR ( dBSRmat *A,
                                    dvector *b,
                                    dvector *x,
                                    ITS_param *itparam,
                                    ILU_param *iluparam,
                                    AMG_param *amgparam,
                                    ivector *neigh,
                                    ivector *order)
{
    // --------------------------------------------------------------
    // Part 1: prepare
    // --------------------------------------------------------------
    // parameters of iterative method
    const INT print_level = itparam->print_level;
    const INT max_levels  = amgparam->max_levels;

    // return variable
    INT status = FASP_SUCCESS;
    
    // local variables
    INT scaled = 1;
   

    // data of AMG for pressure block
    static   AMG_data *mgl;
    //AMG_data *mgl = fasp_amg_data_create(max_levels);
   
    // data for ILU for reservoir block
    static ILU_data LU;
    
    // timing
    REAL setup_start, setup_end, solver_start, solver_end;
    REAL solver_time, setup_time;
    
    // --------------------------------------------------------------
    // Part 2: set up the preconditioner
    // --------------------------------------------------------------
    fasp_gettime(&setup_start);
    
    if (fasp_called_times == 0) mgl = fasp_amg_data_create(max_levels);
    else if (resetup_c == 1){
        fasp_amg_data_free(mgl,amgparam);
        mgl = fasp_amg_data_create(max_levels);
    }


    // AMG setup for Pressure block: initialize A, b, x for mgl[0]
	if ((fasp_called_times==0)||(resetup_c ==1)) {

    mgl[0].A = fasp_dbsr_getPP(A);
    mgl[0].b = fasp_dvec_create(mgl[0].A.row);
    mgl[0].x = fasp_dvec_create(mgl[0].A.col);
    
    switch (amgparam->AMG_type) {
        case UA_AMG: // Unsmoothed Aggregation AMG
            status = fasp_amg_setup_ua(mgl, amgparam);
            break;
        case SA_AMG: // Smoothed Aggregation AMG
            status = fasp_amg_setup_sa(mgl, amgparam);
            break;
        default: // CLASSIC_AMG
            status = fasp_amg_setup_rs(mgl, amgparam);
            break;
    }
    if (status < 0) goto FINISHED;
    }
    
	// ILU decomposition for the whole reservoir block;

#ifdef _OPENMP

#if ILU_MC_OMP
    // printf("ILU_MC_OMP = %d\n", ILU_MC_OMP);
    fasp_ilu_data_free(&LU);
    if ((status = fasp_ilu_dbsr_setup_mc_omp(A, &mgl[0].A, &LU,iluparam))<0) goto FINISHED;
#else
    if ((fasp_called_times == 0) || (resetup_c == 1)) {
        fasp_ilu_data_free(&LU);
        if ( (status = fasp_ilu_dbsr_setup_levsch_step(A,&LU,iluparam,1))<0 ) goto FINISHED;  
    }

    if ((status = fasp_ilu_dbsr_setup_levsch_step(A,&LU,iluparam,2))<0) goto FINISHED;

#endif

#else
    if ((fasp_called_times == 0) || (resetup_c == 1)) {
        fasp_ilu_data_free(&LU);
        if ( (status = fasp_ilu_dbsr_setup_step(A,&LU,iluparam,1))<0 ) goto FINISHED;  
    }

    if ( (status = fasp_ilu_dbsr_setup_step(A,&LU,iluparam,2))<0 ) goto FINISHED;
#endif
    // check iludata
    if ( (status = fasp_mem_iludata_check(&LU))<0 ) goto FINISHED;
    
    fasp_gettime(&setup_end);
    
    setup_time = setup_end - setup_start;
    
    precond_blkoil_bbsr_data pcdata;
    pcdata.maxit           = amgparam->maxit;
    pcdata.tol             = amgparam->tol;
    pcdata.cycle_type      = amgparam->cycle_type;
    pcdata.smoother        = amgparam->smoother;
    pcdata.smooth_order    = amgparam->smooth_order;
    pcdata.presmooth_iter  = amgparam->presmooth_iter;
    pcdata.postsmooth_iter = amgparam->postsmooth_iter;
    pcdata.coarsening_type = amgparam->coarsening_type;
    pcdata.coarse_solver   = amgparam->coarse_solver;
    pcdata.relaxation      = amgparam->relaxation;
    pcdata.coarse_scaling  = amgparam->coarse_scaling;
    pcdata.max_levels      = mgl[0].num_levels;
    pcdata.mgl_data        = mgl;
    pcdata.scaled          = scaled;
    pcdata.RR              = A;
    pcdata.LU              = &LU;
    pcdata.order           = order;
    pcdata.neigh           = neigh;
    
    precond prec; prec.data=&pcdata;
    prec.fct = fasp_precond_dbsr_CPR;
    
    //--------------------------------------------------------------
    // Part 3: solver
    //--------------------------------------------------------------
    fasp_gettime(&solver_start);
    status       = fasp_solver_dbsr_itsolver(A,b,x,&prec,itparam);
    fasp_gettime(&solver_end);
    
    solver_time = solver_end - solver_start;
    total_setup_time += setup_time ; /**< Total setup time                  */
    total_solve_time += solver_time; /**< Total setup time                  */
    
    if ( print_level > PRINT_NONE ) {
        printf("Setup costs %f seconds.\n", setup_time);
        printf("Iterative solver (maxit: %d) %d costs %f seconds : %5.1f\n",itparam->maxit, status, solver_time,setup_time/(solver_time/status) );
        printf("BSR_Krylov_ASCPR totally costs %f seconds.\n", setup_time + solver_time);
    }
    
FINISHED:

    // fasp_amg_data_free(mgl,amgparam);
    // fasp_ilu_data_free(&LU);
    
    if (status == ERROR_ALLOC_MEM) goto MEMORY_ERROR;
    return status;
    
MEMORY_ERROR:
    printf("### ERROR: Setup failed in %s!\n", __FUNCTION__);
    exit(status);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
