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
    const int lda = m;              // leading dimension of A
    const int nb = atoi(argv[2]);   // tile size
    const int bs = atoi(argv[3]);   // inner blocking size

    assert(m >= nb);
    assert(nb >= bs);
	
    #ifdef DEBUG
    cout << "m = " << m << ", nb = " << nb << ", bs = " << bs << endl;
    cout << "# of threads = " << omp_get_max_threads() << endl;
    #endif

    //////////////////////////////////////////////////////////////////////
    // Allocate memory for matrices
    double *A = new double [m * m];  // A: matrix
    Gen_rand_mat(20240417, m, m, A, lda);

   	const int p = (m % nb == 0) ? m/nb : m/nb+1;  // # of tile rows
	
    double *AT[p][p];  // AT: Tiled matrix
    double *TT[p][p];  // TT: Tiled Householder matrix
    for (int j=0; j<p; j++)
    {
        const int jb = min(m-j*nb,nb);     // tile col size

        #pragma omp parallel for
        for (int i=0; i<p; i++)
        {
            const int ib = min(m-i*nb,nb); // tile row size

            AT[i][j] = new double [ib * jb];
            TT[i][j] = new double [bs * jb];

            // Copy elements of A to AT
            const double* Aij = A+(j*nb*lda + i*nb);
            for (int jj=0; jj<jb; jj++)
                for (int ii=0; ii<ib; ii++)
                        AT[i][j][ ii + jj*ib ] = Aij[ ii + jj*lda ];
        }
    }

    // Timer start
    double time = omp_get_wtime();

    //////////////////////////////////////////////////////////////////////
    // Main loop
    #pragma omp parallel
    {
        #pragma omp master
        {
            for (int k=0; k<p; k++)
            {
                const int kb = min(m-k*nb,nb);     // tile row size

                #pragma omp task depend(inout:AT[k][k],TT[k][k])
                GEQRT(kb, kb, bs, 
                    AT[k][k], kb, TT[k][k], bs);

                for (int j=k+1; j<p; j++)
                {
                    const int jb = min(m-j*nb,nb); // tile col size

                    #pragma omp task depend(in:AT[k][k],TT[k][k]) depend(inout:AT[k][j])
                    LARFB(PlasmaLeft, PlasmaTrans, 
                        kb, jb, kb, bs, 
                        AT[k][k], kb, TT[k][k], bs, AT[k][j], kb);
                }

                for (int i=k+1; i<p; i++)
                {
                    const int ib = min(m-i*nb,nb); // tile row size

                    #pragma omp task depend(in:AT[k][k],TT[k][k]) depend(inout:AT[k][k],AT[i][k],TT[i][k])
                    TSQRT(kb, kb, ib, kb, bs, 
                        AT[k][k], kb, AT[i][k], ib, TT[i][k], bs);
                    
                    for (int j=k+1; j<p; j++)
                    {
                        const int jb = min(m-j*nb,nb); // tile col size

                        #pragma omp task depend(in:AT[i][k],TT[i][k]) depend(inout:AT[k][j],AT[i][j])
                        SSRFB(PlasmaLeft, PlasmaTrans, 
                            kb, jb, ib, jb, kb, bs,
                            AT[i][k], ib, TT[i][k], bs, AT[k][j], kb, AT[i][j], ib);
                    }
                }
            }
        }
    }
    // end of main loop
    //////////////////////////////////////////////////////////////////////

    // Timer stop
    time = omp_get_wtime() - time;
    cout << m << ", " << nb << ", " << bs << ", " << time << endl;

    //////////////////////////////////////////////////////////////////////
    #ifdef DEBUG
    double *QT[p][p];          // QT: Tiled matrix

    for (int j=0; j<p; j++)
    {
        const int jb = min(m-j*nb,nb);     // tile col size

        #pragma omp parallel for
        for (int i=0; i<p; i++)
        {
            const int ib = min(m-i*nb,nb); // tile row size
            QT[i][j] = new double [ib * jb];

            if (i==j)
            {
                // Diagonal tiles are set to be identity matrix
                assert(EXIT_SUCCESS == LAPACKE_dlaset(LAPACK_COL_MAJOR, 'g', ib, jb, 0.0, 1.0, QT[i][j], ib));
            }
            else
            {
                // Non-diagonal tiles are set to be zero matrix
                assert(EXIT_SUCCESS == LAPACKE_dlaset(LAPACK_COL_MAJOR, 'g', ib, jb, 0.0, 0.0, QT[i][j], ib));
            }
        }
    }

    // Make otrhogonal tiled matrix QT from AT and TT
    for (int k = p-1; k>=0; k--)
    {
        const int kb = min(m-k*nb,nb);     // tile row size

        for (int i=p-1; i>k; i--)
        {
            const int ib = min(m-i*nb,nb); // tile row size

            for (int j=k; j<p; j++)
            {
                const int jb = min(m-j*nb,nb); // tile col size

                SSRFB(PlasmaLeft, PlasmaNoTrans, 
                    kb, jb, ib, jb, kb, bs,
                    AT[i][k], ib, TT[i][k], bs, QT[k][j], kb, QT[i][j], ib);
            }
        }
        for (int j=k; j<p; j++)
        {
            const int jb = min(m-j*nb,nb); // tile col size

            LARFB(PlasmaLeft, PlasmaNoTrans, 
                kb, jb, kb, bs, 
                AT[k][k], kb, TT[k][k], bs, QT[k][j], kb);
        }
    }

    //////////////////////////////////////////////////////////////////////
    // Check orthogonality of QT
    double *Q = new double [m * m]; // Q: orthogonal matrix
    for (int j=0; j<p; j++)
    {
        const int jb = min(m-j*nb,nb);     // tile col size

        #pragma omp parallel for
        for (int i=0; i<p; i++)
        {
            const int ib = min(m-i*nb,nb); // tile row size

            // Copy elements of QT to Q
            double* Qij = Q+(j*nb*lda + i*nb);
            for (int jj=0; jj<jb; jj++)
                for (int ii=0; ii<ib; ii++)
                        Qij[ ii + jj*lda ] = QT[i][j][ ii + jj*ib ];
        }
    }

    double *Id = new double [m * m]; // I: identity matrix
    assert(EXIT_SUCCESS == LAPACKE_dlaset_work(LAPACK_COL_MAJOR, 'g', m, m, 0.0, 1.0, Id, lda));

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
    double *R = new double [m * m]; // R: upper triangular matrix

    // Copy the upper triangular tiles of AT to R
    for (int j=0; j<p; j++)
    {
        const int jb = min(m-j*nb,nb);     // tile col size

        #pragma omp parallel for
        for (int i=0; i<p; i++)
        {
            const int ib = min(m-i*nb,nb); // tile row size
            double* Rij = R+(j*nb*lda + i*nb);

            if (j>=i)          // lower part
            {
                if (j==i)
                {
                    for (int jj=0; jj<jb; jj++)
                        for (int ii=0; ii<ib; ii++)
                            if (jj >= ii)
                                Rij[ ii + jj*lda ] = AT[i][j][ ii + jj*ib ];
                            else
                                Rij[ ii + jj*lda ] = 0.0;
                }
                else
                {
                    for (int jj=0; jj<jb; jj++)
                        for (int ii=0; ii<ib; ii++)
                            Rij[ ii + jj*lda ] = AT[i][j][ ii + jj*ib ];
                }
            }
            else
            {
                for (int jj=0; jj<jb; jj++)
                        for (int ii=0; ii<ib; ii++)
                            Rij[ ii + jj*lda ] = 0.0;
            }
        }
    }

    for (int j=0; j<p; j++)
    {
        const int jb = min(m-j*nb,nb);     // tile col size

        #pragma omp parallel for
        for (int i=0; i<p; i++)
        {
            const int ib = min(m-i*nb,nb); // tile row size

            // Copy elements of QT to Q
            double* Qij = Q+(j*nb*lda + i*nb);
            for (int jj=0; jj<jb; jj++)
                for (int ii=0; ii<ib; ii++)
                        Qij[ ii + jj*lda ] = QT[i][j][ ii + jj*ib ];
        }
    }

    // |A|_oo
    double normA = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', m, m, A, lda, work);

    // Compute A - Q*R.
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, m, m, m, 1.0, Q, lda, R, lda, -1.0, A, lda);

    // |A - Q*R|_oo
    double error = LAPACKE_dlange_work(LAPACK_COL_MAJOR, 'I', m, m, A, lda, work);

    // normalize the result
    // |A-QR|_oo / (|A|_oo * m)
    // error /= (normA * m);

    cout << "Accuracy: " << error << endl;

    delete [] work;
    delete [] Q;

    for (int j=0; j<p; j++)
        for (int i=0; i<p; i++)
        {
            delete [] QT[i][j];
        }

    #endif
    //////////////////////////////////////////////////////////////////////

    for (int j=0; j<p; j++)
        for (int i=0; i<p; i++)
        {
            delete [] AT[i][j];
            delete [] TT[i][j];
        }
    delete [] A;

    return EXIT_SUCCESS;
}
