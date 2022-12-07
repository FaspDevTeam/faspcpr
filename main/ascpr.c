/*! \file  main.c
 *
 *  \brief An Improved CPR-type Preconditioner Test Function for the Black Oil Model
 * 
 *---------------------------------------------------------------------------------
 *  Copyright (C) 2015--2022 by the FASP team. All rights reserved.
 *---------------------------------------------------------------------------------
 */

#include <time.h>
#include <stdlib.h>

#include "fasp.h"
#include "fasp_functs.h"
#include "faspcpr.h"
#include "faspcpr_functs.h"


/**
 * \fn int main (int argc, const char * argv[])
 *
 * \brief This is the main function for an improved CPR-type preconditioner test.
 *
 * \author Li Zhao, Chunsheng Feng, Chensong Zhang
 * \date   01/08/2022
 *
 */
int main(int argc, const char *argv[])
{
    const INT  print_level = 3;    // how much information to print out (0~10)
    const REAL tolerance   = 1e-5; // tolerance for accepting the solution
    const INT  mu = 0;             // given a threshold (mu >= 0),  the mu determines whether a new "setup phase" is necessary,
    // Parameters
    ITS_param  itparam;  // parameters for itsolver
    AMG_param  amgparam; // parameters for AMG
    ILU_param  iluparam; // parameters for ILU
    
    // Matrix (BSR format)
    dBSRmat Absr;

    // Right-hand side and solutions
    dvector    b, x;
    
    // Local variables
    int i, N = 1;
    int problem_num;  // test problem number
    time_t lt  = time(NULL);
    

    // Initialize input parameters
    fasp_param_amg_init(&amgparam);
    fasp_param_solver_init(&itparam);
    fasp_param_ilu_init(&iluparam);
    

    // parameters for itsolver
    fasp_param_solver_init(&itparam);
    itparam.print_level   = print_level;
    itparam.tol           = tolerance;
    itparam.maxit         = 100;
    itparam.itsolver_type = SOLVER_VGMRES;
    itparam.restart       = 28; // the number of restart for GMRES


    // parameters for ILU
    iluparam.ILU_lfil             = 0; // fill-level for ILU(k), e.g., ILU(0) <==> ILU_lfil = 0.
    iluparam.print_level          = print_level;
	

    // parameters for AMG
    amgparam.print_level          = print_level;
    amgparam.maxit                = 1;
    amgparam.coarse_dof           = 10000;
    // amgparam.smooth_order         = NO_ORDER; 
    amgparam.AMG_type             = UA_AMG;
    amgparam.cycle_type           = NL_AMLI_CYCLE;
    amgparam.aggregation_type     = NPAIR; // non-symmetric aggregation!

    // We should choose one of these four direct solvers for coarse spaces
#if WITH_UMFPACK
    if (!fasp_called_times) printf("\namgparam.coarse_solver = umfpack\n"); 
    amgparam.coarse_solver        = SOLVER_UMFPACK;
#endif
        
#if WITH_SUPERLU
    if (!fasp_called_times) printf("\namgparam.coarse_solver = superlu\n");
    amgparam.coarse_solver        = SOLVER_SUPERLU;
#endif
        
#if WITH_MUMPS
    if (!fasp_called_times) printf("\namgparam.coarse_solver = mumps\n"); 
    amgparam.coarse_solver        = SOLVER_MUMPS;
#endif

#if WITH_PARDISO
    if (!fasp_called_times) printf("\namgparam.coarse_solver = pardiso\n");
    amgparam.coarse_solver        = SOLVER_PARDISO;
#endif


    printf("\n\n");
    printf("------------------------- Test starts at -------------------------\n");
    printf("%s",asctime(localtime(&lt))); // output starting local time
    printf("------------------------------------------------------------------\n");
    
    
    // Problem 1. Black oil in BSR format
    problem_num = 1;

    printf("\n=====================================================\n");
    printf("Test Problem Number %d ...\n", problem_num);
    printf("1200X2200X170 Three-Phase Black Oil in BSR format");
    printf("\n=====================================================\n");
    
    if (problem_num == 1) {  
        char *matrixfile1 = "data/A_test1.dat";
        char *rhsfile1    = "data/b_test1.dat";

		// Read matrix from file "matrixfile1"
		fasp_dbsr_read(matrixfile1, &Absr);
		// Read right hand side from file "rhsfile1"
		fasp_dvec_read (rhsfile1, &b);
	}
	
    // Allocate space for initial solution
    fasp_dvec_alloc(b.row, &x);
        
    // ASCPR method
    printf("------------------------------------------------------------------\n");
    // printf("For a linear algebraic system (BSR format) arising from the black oil model \n");
    printf("An Adaptive SETUP CPR preconditioner (ASCPR) is tested, N = %d ...\n", N);

    // start time
#ifdef _OPENMP
    REAL omp_start_time = omp_get_wtime();
#else
    clock_t start_time = clock();
#endif

    for(i = 0; i < N; i++){
	printf("The number of ASCPR-GMRES calls: %d\n", i);
       	// reset initial guess
       	fasp_dvec_set(b.row, &x, 0.0); 
       	
       	// call FASP_BSRSOL_ASCPR, i.e., ASCPR-GMRES solver
	FASP_BSRSOL_ASCPR(&Absr, &b, &x, &itparam, &iluparam, &amgparam, mu);
    }

    // end time
#ifdef _OPENMP
    REAL omp_end_time = omp_get_wtime();
#else
    clock_t end_time = clock();
#endif

    // print time
#ifdef _OPENMP
    printf("ASCPR-GMRES-OMP time: %.4f\n", omp_end_time - omp_start_time);
#else
    printf("ASCPR-GMRES-SEQ time: %.4f\n", (REAL)(end_time - start_time) / (REAL)CLOCKS_PER_SEC);
#endif
    // Clean up memory
    fasp_dvec_free(&x);
    fasp_dbsr_free(&Absr);
    fasp_dvec_free(&b);

  
    return FASP_SUCCESS;
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
