//
//  main.cpp
//
//  Created by T. Suzuki on 2024/04/17.
//

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <omp.h>
#include <mkl.h>
#include <plasma_core_blas.h>

#include "CoreBlas.h"

using namespace std;

// Generate random matrix
void Gen_rand_mat(const long int seed, const int m, const int n, double *A, const int lda)
{
	srand(seed);

    for (int j=0; j<n; j++)
	    for (int i=0; i<m; i++)
            A[i+j*lda] = 1.0 - 2*(double)rand() / RAND_MAX;
}

// Show matrix
void Show_mat(const int m, const int n, double *A, const int lda)
{
  	for (int i=0; i<m; i++)
    {
		for (int j=0; j<n; j++)
			printf("% 5.4lf, ", A[i + j*lda]);
		cout << endl;
	}
	cout << endl;
}

// Debug mode
#define DEBUG

// Trace mode
//#define TRACE

#ifdef TRACE
extern void trace_cpu_start();
extern void trace_cpu_stop(const char *color);
extern void trace_label(const char *color, const char *label);
#endif

int main(int argc, const char ** argv)
{

    //Usage "a.out [size of matrix: M ] [tile size: NB] [inner block size: IB]"
    assert(argc > 3);

    const int m =  atoi(argv[1]);   // matrix size
    const int n = m;                // A should be square
    const int lda = m;              // leading dimension of A
    const int mb = atoi(argv[2]);   // tile size
    const int nb = mb;              // square tile
    const int ib = atoi(argv[3]);   // inner blocking size

    assert(m >= mb);
    assert(mb >= ib);
	
    #ifdef DEBUG
    cout << "m = " << m << ", mb = " << mb << ", ib = " << ib << endl;
    cout << "# of threads = " << omp_get_max_threads() << endl;
    #endif

    //////////////////////////////////////////////////////////////////////
    // Allocate memory for matrices
    double *A = new double [m * n]; // A: matrix
    Gen_rand_mat(20240417, m, n, A, lda);

   	const int nmb = (m % mb == 0) ? m/mb : m/mb+1;  // # of tile rows
	const int nnb = (n % nb == 0) ? n/nb : n/nb+1;  // # of tile cols
	
    double *AT[nmb][nnb];          // AT: Tiled matrix
    double *TT[nmb][nnb];          // TT: Tiled Householder matrix
    for (int j=0; j<nnb; j++)
    {
        const int nbj = min(n-j*nb,nb);     // tile col size
        for (int i=0; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size
            AT[i][j] = new double [mbi * nbj];
            TT[i][j] = new double [ib * nbj];

            // Copy elements of A to AT
            const double* Aij = A+(j*nb*lda + i*mb);
            for (int jj=0; jj<nbj; jj++)
                for (int ii=0; ii<mbi; ii++)
                        AT[i][j][ ii + jj*mbi ] = Aij[ ii + jj*lda ];
        }
    }

    // Timer start
    double time = omp_get_wtime();

    //////////////////////////////////////////////////////////////////////
    // Main loop
    for (int k=0; k<min(nmb,nnb); k++)
    {
        const int mbk = min(m-k*mb,mb);     // tile row size

        GEQRT(mbk, mbk, ib, 
            AT[k][k], mbk, TT[k][k], ib);
        // cout << "GEQRT[" << k << "][" << k << "]" << endl;

        for (int j=k+1; j<nnb; j++)
        {
            const int nbj = min(n-j*nb,nb); // tile col size

            LARFB(PlasmaLeft, PlasmaTrans, 
                mbk, nbj, mbk, ib, 
                AT[k][k], mbk, TT[k][k], ib, AT[k][j], mbk);
            // cout << "LARFB[" << k << "][" << j << "]" << endl;
        }

        for (int i=k+1; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size

            TSQRT(mbk, mbk, mbi, mbk, ib, 
                AT[k][k], mbk, AT[i][k], mbi, TT[i][k], ib);
            // cout << "TSQRT[" << i << "][" << k << "]" << endl;
            
            for (int j=k+1; j<nnb; j++)
            {
                const int nbj = min(n-j*nb,nb); // tile col size

                SSRFB(PlasmaLeft, PlasmaTrans, 
                    mbk, nbj, mbi, nbj, mbk, ib,
                    AT[i][k], mbi, TT[i][k], ib, AT[k][j], mbk, AT[i][j], mbi);
                // cout << "SSRFB[" << i << "][" << j << "]" << endl;
            }
        }
    }
    // end of main loop
    //////////////////////////////////////////////////////////////////////

    // Timer stop
    time = omp_get_wtime() - time;
    cout << m << ", " << nb << ", " << ib << ", " << time << endl;

    // Show_mat(m, n, A, lda);

    //////////////////////////////////////////////////////////////////////
    #ifdef DEBUG
    double *QT[nmb][nnb];          // QT: Tiled matrix

    // cout << "nmb = " << nmb << ", nnb = " << nnb << endl;
    for (int j=0; j<nnb; j++)
    {
        const int nbj = min(n-j*nb,nb);     // tile col size
        for (int i=0; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size
            QT[i][j] = new double [mbi * nbj];

            if (i==j)
            {
                // Diagonal tiles are set to be identity matrix
                LAPACKE_dlaset(LAPACK_COL_MAJOR, 'g', mbi, nbj, 0.0, 1.0, QT[i][j], mbi);
            }
            else
            {
                // Non-diagonal tiles are set to be zero matrix
                LAPACKE_dlaset(LAPACK_COL_MAJOR, 'g', mbi, nbj, 0.0, 0.0, QT[i][j], mbi);
            }
        }
    }

    // Make otrhogonal tiled matrix QT from AT and TT
    for (int k = min(nmb,nnb)-1; k>=0; k--)
    {
        const int mbk = min(m-k*mb,mb);     // tile row size

        for (int i=nmb-1; i>k; i--)
        {
            const int mbi = min(m-i*mb,mb); // tile row size

            for (int j=k; j<nnb; j++)
            {
                const int nbj = min(n-j*nb,nb); // tile col size

                SSRFB(PlasmaLeft, PlasmaNoTrans, 
                    mbk, nbj, mbi, nbj, mbk, ib,
                    AT[i][k], mbi, TT[i][k], ib, QT[k][j], mbk, QT[i][j], mbi);
                // cout << "SSRFB[" << i << "][" << j << "]" << endl;
            }
        }
        for (int j=k; j<nnb; j++)
        {
            const int nbj = min(n-j*nb,nb); // tile col size

            LARFB(PlasmaLeft, PlasmaNoTrans, 
                mbk, nbj, mbk, ib, 
                AT[k][k], mbk, TT[k][k], ib, QT[k][j], mbk);
            // cout << "LARFB[" << k << "][" << j << "]" << endl;
        }
    }

    //////////////////////////////////////////////////////////////////////
    // Check orthogonality of QT
    double *Q = new double [m * m]; // Q: orthogonal matrix
    for (int j=0; j<nnb; j++)
    {
        const int nbj = min(n-j*nb,nb);     // tile col size
        for (int i=0; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size

            // Copy elements of QT to Q
            double* Qij = Q+(j*nb*lda + i*mb);
            for (int jj=0; jj<nbj; jj++)
                for (int ii=0; ii<mbi; ii++)
                        Qij[ ii + jj*lda ] = QT[i][j][ ii + jj*mbi ];
        }
    }

    double *Id = new double [m * m]; // I: identity matrix
    LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'g', m, m, 0.0, 1.0, Id, lda);

    // Perform Id - Q^T * Q
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, m, m, 1.0, Q, lda, -1.0, Id, lda);

    // work array of size m is needed for computing L_oo norm
    double *work = new double [m];

    // |I - Q^T * Q|_oo
    double ortho = LAPACKE_dlansy_work(LAPACK_COL_MAJOR, 'I', 'u', m, Id, m, work);

    // normalize the result
    // |Id - Q^T * Q|_oo / n
    // ortho /= m;
    cout << "Orthogonality: " << ortho << endl;

    delete [] Id;

    //////////////////////////////////////////////////////////////////////
    // Check the accuracy of A - Q * R

    // Upper triangular matrix R.
    double *R = new double [m * n]; // R: upper triangular matrix

    // Copy the upper triangular tiles of AT to R
    for (int j=0; j<nnb; j++)
    {
        const int nbj = min(n-j*nb,nb);     // tile col size

        for (int i=0; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size
            double* Rij = R+(j*nb*lda + i*mb);

            if (j>=i)          // lower part
            {
                if (j==i)
                {
                    // cout << "(" << i << "," << j << ")" << endl;
                    for (int jj=0; jj<nbj; jj++)
                        for (int ii=0; ii<mbi; ii++)
                            if (jj >= ii)
                                Rij[ ii + jj*lda ] = AT[i][j][ ii + jj*mbi ];
                            else
                                Rij[ ii + jj*lda ] = 0.0;
                }
                else
                {
                    // cout << "(" << i << "," << j << ")" << endl;
                    for (int jj=0; jj<nbj; jj++)
                        for (int ii=0; ii<mbi; ii++)
                            Rij[ ii + jj*lda ] = AT[i][j][ ii + jj*mbi ];
                }
            }
            else
            {
                // cout << "(" << i << "," << j << ")" << endl;
                for (int jj=0; jj<nbj; jj++)
                        for (int ii=0; ii<mbi; ii++)
                            Rij[ ii + jj*lda ] = 0.0;
            }
        }
    }

    for (int j=0; j<nnb; j++)
    {
        const int nbj = min(n-j*nb,nb);     // tile col size
        for (int i=0; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size

            // Copy elements of QT to Q
            double* Qij = Q+(j*nb*lda + i*mb);
            for (int jj=0; jj<nbj; jj++)
                for (int ii=0; ii<mbi; ii++)
                        Qij[ ii + jj*lda ] = QT[i][j][ ii + jj*mbi ];
        }
    }

    // |A|_oo
    double normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', m, n, A, lda, work);

    // Compute A - Q*R.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, n, m, 1.0, Q, lda, R, lda, -1.0, A, lda);

    // |A - Q*R|_oo
    double error = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', m, n, A, lda, work);

    // normalize the result
    // |A-QR|_oo / (|A|_oo * n)
    // error /= (normA * n);

    cout << "Accuracy: " << error << endl;

    delete [] work;
    delete [] Q;

    for (int j=0; j<nnb; j++)
        for (int i=0; i<nmb; i++)
        {
            delete [] QT[i][j];
        }

    #endif
    //////////////////////////////////////////////////////////////////////

    for (int j=0; j<nnb; j++)
        for (int i=0; i<nmb; i++)
        {
            delete [] AT[i][j];
            delete [] TT[i][j];
        }
    delete [] A;

    return EXIT_SUCCESS;
}
