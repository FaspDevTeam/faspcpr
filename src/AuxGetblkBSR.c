/*! \file  AuxGetblkBSR.c
 *
 *  \brief Get sub-blocks of dBSRmat matrices
 *
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

/*---------------------------------*/
/*--      Public Functions       --*/
/*---------------------------------*/

/**
 * \fn dCSRmat fasp_dbsr_getPP (dBSRmat *A)
 *
 * \brief get pressure block from the reservoir block
 *
 * \param A    pointer to the BSR format matrix
 *
 * \return     dCSRmat Pressure matrix if succeed, NULL if fail
 *
 * \author Xiaozhe Hu
 * \date 12/25/2010
 
 * \note Works for general nb (Xiaozhe)
 *
 * Modified by Chunsheng Feng, Zheng Li on 08/28/2012
 */
dCSRmat fasp_dbsr_getPP (dBSRmat *A)
{
	// information about A
	const INT ROW = A->ROW;
	const INT COL = A->COL;
	const INT NNZ = A->NNZ;
	const INT nc = A->nb;
	const INT nc2 = nc*nc;
	
	REAL *val = A->val;
	INT *IA = A->IA;
	INT *JA = A->JA;
	
	// Pressure block
	dCSRmat P_csr = fasp_dcsr_create(ROW, COL, NNZ);
	REAL *Pval=P_csr.val;
	
	// local variable
	INT i,j;
	INT status = FASP_SUCCESS;
    
#ifdef _OPENMP
    // variables for OpenMP
    INT myid, mybegin, myend;
    INT nthreads = fasp_get_num_threads();
#endif
    
	// get pressure block
	memcpy(P_csr.JA, JA, NNZ*sizeof(INT));
	memcpy(P_csr.IA, IA, (ROW+1)*sizeof(INT));
    
#ifdef _OPENMP
	if (NNZ > OPENMP_HOLDS) {
#pragma omp parallel for private (myid, mybegin, myend, i)
	    for ( myid = 0; myid < nthreads; myid++ ) {
            fasp_get_start_end(myid, nthreads, NNZ, &mybegin, &myend);
            for ( i = mybegin; i < myend; i++ ) {
                Pval[i] = val[i*nc2];
            }
        }
	}
	else {
#endif
	    //for ( i = NNZ, j = NNZ*nc2-nc2; i--; j-=nc2 ) {
        for ( i = NNZ, j = NNZ*nc2-nc2; i--; j-=nc2 ) {
            Pval[i] = val[j];
	    }
        
#ifdef _OPENMP
    }
#endif
    
	// compress CSR format
	status = fasp_dcsr_compress_inplace(&P_csr,1e-8);
	if ( status < 0 ) goto FINISH_ERROR;
	
	// return P
	return P_csr;
	
FINISH_ERROR:
	printf("### ERROR: Cannot compress the matrix P %s!\n", __FUNCTION__);
	exit(ERROR_ALLOC_MEM);
}

/*---------------------------------*/
/*--        End of File          --*/
/*---------------------------------*/
