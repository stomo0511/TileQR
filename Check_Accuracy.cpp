//
//  Check_Accuracy
//
//  Created by T. Suzuki on 2013/08/16.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <iostream>
#include "Check_Accuracy.h"

using namespace std;

void Check_Accuracy( const int M, const int N, double *mA, double *mQ, double *mR )
{
    //////////////////////////////////////////////////////////////////////
    // Check Orthogonarity

    // Set Id to the identity matrix
    int mn = min(M,N);
    double* Id = new double[ mn * mn ];
    for (int i=0; i<mn; i++)
        for (int j=0; j<mn; j++)
        Id[ i + j*mn ] = (i == j) ? 1.0 : 0.0;

    double alpha = 1.0;
    double beta  = -1.0;
        
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, 
            N, M, alpha, mQ, M, beta, Id, N);
        
    double* Work = new double[ mn ];
    double ortho = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, 'I', 'U', 
                        mn, Id, mn, Work);
    delete [] Work;

    // normalize the result
    // |Id - Q^T * Q|_oo / n
    // ortho /= mn;
    std::cout << "norm(I-Q'*Q) = " << ortho << std::endl;

    // Check Orthogonarity END
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Check Residure
    Work = new double[ M ];

    // |A|_oo
    double normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', M, N, mA, M, Work);

    alpha = -1.0;
    beta  =  1.0;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
            M, N, M, alpha, mQ, M, mR, M, beta, mA, M);
        
    double normQ = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', 
                    M, N, mA, M, Work);

    // normalize the result
    // |A-QR|_oo / (|A|_oo * n)
    // normQ /= (normA * N);
    std::cout << "norm(A-Q*R) = " << normQ << std::endl;

    // Check Residure END
    //////////////////////////////////////////////////////////////////////

    delete [] Id;
    delete [] Work;
    // delete [] QR;
}

