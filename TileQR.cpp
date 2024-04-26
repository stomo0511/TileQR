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

// Generate random matrix
//
// @param[in] seed: random seed
// @param[in] m: row size of matrix
// @param[in] n: col size of matrix
// @param[out] A: (m x n) matrix
// @param[in] lda: leading dimension of A
//
void Gen_rand_mat(const long int seed, const int m, const int n, double *A, const int lda)
{
	srand(seed);

    for (int j=0; j<n; j++)
	    for (int i=0; i<m; i++)
            A[i+j*lda] = 1.0 - 2*(double)rand() / RAND_MAX;
}

// Show matrix
//
// @param[in] m: row size of matrix
// @param[in] n: col size of matrix
// @param[in] A: (m x n) matrix
// @param[in] lda: leading dimension of A
//
void Show_mat(const int m, const int n, double *A, const int lda)
{
  	for (int i=0; i<m; i++)
    {
		for (int j=0; j<n; j++)
			printf("% 5.4lf, ", A[i + j*lda]);
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

// Convert column major format to Tiled layout
//
void cm2ccrb(const int m, const int n, const int lda, const double *A, const int mb, const int nb, double* B)
{
   	const int nmb = (m % mb == 0) ? m/mb : m/mb+1;  // # of tile rows
	const int nnb = (n % nb == 0) ? n/nb : n/nb+1;  // # of tile cols

    for (int j=0; j<nnb; j++)
    {
        int jb = std::min(n-j*nb,nb);
        for (int i=j; i<nmb; i++)
        {
            int ib = std::min(m-i*mb,mb);
            const double* Aij = A+(j*nb*lda + i*mb);
            double* Bij = B+(j*nb*m + i*mb*jb);

            #pragma omp task depend(in: Aij[0:m*jb]) depend(out: Bij[0:ib*jb])
            {
                #ifdef TRACE
                trace_cpu_start();
                trace_label("Yellow", "Conv.");
                #endif

                for (int jj=0; jj<jb; jj++)
                    for (int ii=0; ii<ib; ii++)
                        Bij[ ii + jj*ib ] = Aij[ ii + jj*lda ];

                #ifdef TRACE
                trace_cpu_stop("Yellow");
                #endif
            }
        }
    }
}

void ccrb2cm(const int m, const int n, const int lda, double *A, const int mb, const int nb, const double* B)
{
   	const int nmb = (m % mb == 0) ? m/mb : m/mb+1;  // # of tile rows
	const int nnb = (n % nb == 0) ? n/nb : n/nb+1;  // # of tile cols

    for (int j=0; j<nnb; j++)
    {
        int jb = std::min(m-j*nb,nb);
        for (int i=j; i<nmb; i++)
        {
            int ib = std::min(m-i*mb,mb);
            double* Aij = A+(j*nb*lda + i*nb);
            const double* Bij = B+(j*nb*lda + i*nb*jb);

            #pragma omp task depend(in: Bij[0:ib*jb]) depend(out: Aij[0:m*jb])
            {
                #ifdef TRACE
                trace_cpu_start();
                trace_label("Violet", "Conv.");
                #endif

                for (int jj=0; jj<jb; jj++)
                    for (int ii=0; ii<ib; ii++)
                        Aij[ ii+jj*lda ] = Bij[ ii+jj*ib ];

                #ifdef TRACE
                trace_cpu_stop("Violet");
                #endif
            }
        }
    }
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

    double *A = new double [m * n]; // A: matrix
    Gen_rand_mat(20240417, m, n, A, lda);

   	const int nmb = (m % mb == 0) ? m/mb : m/mb+1;  // # of tile rows
	const int nnb = (n % nb == 0) ? n/nb : n/nb+1;  // # of tile cols

    assert(m >= mb);
    assert(mb >= ib);
	
    #ifdef DEBUG
    std::cout << "m = " << m << ", mb = " << mb << ", ib = " << ib << std::endl;
    std::cout << "# of threads = " << omp_get_max_threads() << std::endl;
    #endif
	
    double *AT[nmb][nnb];          // AT: Tiled matrix
    double *TT[nmb][nnb];          // TT: Tiled Householder matrix
    for (int j=0; j<nnb; j++)
    {
        const int nbj = std::min(n-j*nb,nb);     // tile col size
        for (int i=0; i<nmb; i++)
        {
            const int mbi = std::min(m-i*mb,mb); // tile row size
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

    for (int k=0; k<std::min(nmb,nnb); k++)
    {
        const int mbk = std::min(m-k*mb,mb);     // tile row size

        GEQRT(mbk, mbk, ib, 
            AT[k][k], mbk, TT[k][k], ib);
        std::cout << "GEQRT[" << k << "][" << k << "]" << std::endl;

        for (int j=k+1; j<nnb; j++)
        {
            const int nbj = std::min(n-j*nb,nb); // tile col size

            LARFB(PlasmaLeft, PlasmaTrans, 
                mbk, nbj, mbk, ib, 
                AT[k][k], mbk, TT[k][k], ib, AT[k][j], mbk);
            std::cout << "LARFB[" << k << "][" << j << "]" << std::endl;
        }

        for (int i=k+1; i<nmb; i++)
        {
            const int mbi = std::min(m-i*mb,mb); // tile row size

            TSQRT(mbk, mbk, mbi, mbk, ib, 
                AT[k][k], mbk, AT[i][k], mbi, TT[i][k], ib);
            std::cout << "TSQRT[" << i << "][" << k << "]" << std::endl;
            
            for (int j=k+1; j<nnb; j++)
            {
                const int nbj = std::min(n-j*nb,nb); // tile col size

                SSRFB(PlasmaLeft, PlasmaTrans, 
                    mbk, nbj, mbi, nbj, mbk, ib,
                    AT[i][k], mbi, TT[i][k], ib, AT[k][j], mbk, AT[i][j], mbi);
                std::cout << "SSRFB[" << i << "][" << j << "]" << std::endl;
            }
        }
    }

    // Timer stop
    time = omp_get_wtime() - time;
    std::cout << m << ", " << nb << ", " << ib << ", " << time << std::endl;

    for (int j=0; j<nnb; j++)
    {
        const int nbj = std::min(n-j*nb,nb);     // tile col size
        for (int i=0; i<nmb; i++)
        {
            const int mbi = std::min(m-i*mb,mb); // tile row size

            // Copy elements of AT to A
            double* Aij = A+(j*nb*lda + i*mb);
            for (int jj=0; jj<nbj; jj++)
                for (int ii=0; ii<mbi; ii++)
                        Aij[ ii + jj*lda ] = AT[i][j][ ii + jj*mbi ];
        }
    }

    Show_mat(m, n, A, lda);

    #ifdef DEBUG
    #endif

    for (int j=0; j<nnb; j++)
        for (int i=0; i<nmb; i++)
        {
            delete [] AT[i][j];
            delete [] TT[i][j];
        }
    delete [] A;

    return EXIT_SUCCESS;
}
