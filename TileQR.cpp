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
#include <mpi.h>

#include "CoreBlas.h"
#include "mpi-detach.h"

using namespace std;

// Set global variables
int myrank;  // myrank
int nprocs;  // number of MPI process
int nthred;  // number of threads

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

// Correspondence between block and rank
void get_block_rank(int *block_rank, int nt)
{
    int row, col;
    row = col = nprocs;
    int* blocks = new int [nprocs];

    for (int i=0; i<nprocs; i++)
        blocks[i]=0;

    if (nprocs != 1)
	{
        while (1)
		{
            row = row / 2;
            if (row * col == nprocs)
				break;
            col = col / 2;
            if (row * col == nprocs)
				break;
			assert(row*col >= nprocs);
        }
    }

    int tmp_rank = 0;
	int offset = 0;
    for (int j = 0; j < nt; j++)
	{
        for (int i = 0; i < nt; i++)
		{
            block_rank[i + j*nt] = tmp_rank + offset;
            blocks[tmp_rank + offset]++;
            tmp_rank++;
            if (tmp_rank >= col)
				tmp_rank = 0;
        }
        tmp_rank = 0;
        offset = (offset + col >= nprocs) ? 0 : offset + col;
    }

	delete [] blocks;
}

#ifdef TRACE
extern void trace_cpu_start();
extern void trace_cpu_stop(const char *color);
extern void trace_label(const char *color, const char *label);
#endif

int main(int argc, char ** argv)
{
   	// MPI Initialize
	int provided;
	MPI_Init_thread(&argc, &argv, MPI_THREAD_MULTIPLE, &provided);
	assert (MPI_THREAD_MULTIPLE == provided);

    MPI_Comm_rank(MPI_COMM_WORLD, &myrank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
	assert(nprocs % 2 == 0);

	#ifdef TRACE
	trace_getrank(myrank);
	#endif

    //Usage "a.out [size of matrix: M ] [tile size: NB] [inner block size: BS]"
    assert(argc > 3);

    const int m =  atoi(argv[1]);   // matrix size
    const int lda = m;              // leading dimension of A
    const int nb = atoi(argv[2]);   // tile size
    const int bs = atoi(argv[3]);   // inner blocking size

    assert(m >= nb);
    assert(nb >= bs);
	
    #ifdef DEBUG
    if (myrank == 0)
    {
        cout << "m = " << m << ", nb = " << nb << ", bs = " << bs << endl;
        cout << "# of threads = " << omp_get_max_threads() << endl;
    }
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

    // Set block rank
	int* block_rank = new int [p*p];
	get_block_rank(block_rank, p);

    // Timer start
    double time;
    if (myrank == 0)
        time = omp_get_wtime();

    //////////////////////////////////////////////////////////////////////
    // Main loop
    #pragma omp parallel
    {
        #pragma omp master
        {
            for (int k=0; k<p; k++)
            {
                const int kb = min(m-k*nb,nb);     // tile row size

                ///////////////////////////////////////////////////////////////
				// GEQRT
                if (block_rank[k + k*p] == myrank)
                {
                    #pragma omp task depend(inout:AT[k][k]) depend(out:TT[k][k])
                    GEQRT(kb, kb, bs, 
                        AT[k][k], kb, TT[k][k], bs);
                    
                    cout << "GEQRT(" << k << "," << k << ") on rank=" << myrank << endl;
                }

                ///////////////////////////////////////////////////////////////
				// GEQRT Send task
                if (block_rank[k + k*p] == myrank)
                {
                    int *send_flags = new int[nprocs];

                    // Send tiles A[k][k] and T[k][k] to j-direction
                    // Initialize send_flags
                    for (int ii=0; ii<nprocs; ii++)
                        send_flags[ii] = 0;
                        
                    for (int jj = k + 1; jj < p; jj++)
                    {
                        if (!send_flags[block_rank[k + jj*p]])
                            send_flags[block_rank[k + jj*p]] = 1;
                    }

                    for (int dst=0; dst<nprocs; dst++)
                    {
                        if (send_flags[dst] && dst != myrank)
                        {
                            #pragma omp task depend(in: AT[k][k], TT[k][k])
                            {
                                double *ATT = new double[(kb+bs)*kb];    // private 宣言で共有変数にできるかも（要確認）
                                
                                cblas_dcopy(kb*kb, AT[k][k], 1, ATT, 1);
                                cblas_dcopy(bs*kb, TT[k][k], 1, ATT+kb*kb, 1);

                                // cout << "Send ATT from " << myrank << " to " << dst << " with tag " << k+k*p << endl;
                                MPI_Send(ATT, (kb+bs)*kb, MPI_DOUBLE, dst, k+k*p, MPI_COMM_WORLD);

                                delete [] ATT;
                            }
                        }
                    }

                    // Send tiles A[k][k] to i-direction
                    // Initialize send_flags
                    for (int ii=0; ii<nprocs; ii++)
                        send_flags[ii] = 0;
                        
                    for (int ii = k + 1; ii < p; ii++)
                    {
                        if (!send_flags[block_rank[ii + k*p]])
                            send_flags[block_rank[ii + k*p]] = 1;
                    }

                    for (int dst=0; dst<nprocs; dst++)
                    {
                        if (send_flags[dst] && dst != myrank)
                        {
                            #pragma omp task depend(in: AT[k][k])
                            {
                                // cout << "Send AT from " << myrank << " to " << dst << " with tag " << 100000+k+k*p << endl;
                                MPI_Send(AT[k][k], kb*kb, MPI_DOUBLE, dst, 100000+(k+k*p), MPI_COMM_WORLD);
                            }
                        }
                    }

                    delete [] send_flags;
                } // GEQRT Send task End

                ///////////////////////////////////////////////////////////////
				// GEQRT Recv task
                if (block_rank[k + k*p] != myrank)
                {
                    int recv_flag = 0;

                    // Recv tiles A[k][k] and T[k][k] from the ranks of tile A[k][k]
					for (int jj = k + 1; jj < p; jj++)
					{
						if (myrank == block_rank[k + jj*p])
							recv_flag = 1;
					}

                    if (recv_flag)
					{
						#pragma omp task depend(out: AT[k][k], TT[k][k])
						{
                            MPI_Status stat;
                            double *ATT = new double[(kb+bs)*kb];    // private 宣言で共有変数にできるかも（要確認）

                            // cout << "Recv ATT from " << block_rank[k + k*p] << " to " << myrank << " with tag " << k+k*p << endl;
							MPI_Recv(ATT, (kb+bs)*kb, MPI_DOUBLE, block_rank[k + k*p], k+k*p, MPI_COMM_WORLD, &stat);
                            
                            cblas_dcopy(kb*kb, ATT, 1, AT[k][k], 1);
                            cblas_dcopy(bs*kb, ATT+kb*kb, 1, TT[k][k], 1);

                            delete [] ATT;
                            // cout << "Recv A[" << k << "][" << k << "], and T[" << k << "," << k << "] on rank=" << myrank << endl;
						}
                    }

   					recv_flag = 0;

                    // Recv tiles A[k][k] from the ranks of tile A[k][k]
 					for (int ii = k + 1; ii < p; ii++)
    				{
	    				if (myrank == block_rank[ii + k*p])
		    				recv_flag = 1;
			    	}

                    if (recv_flag)
					{
						#pragma omp task depend(out: AT[k][k])
						{
							MPI_Status stat;

                            // cout << "Recv AT from " << block_rank[k + k*p] << " to " << myrank << " with tag " << 100000+k+k*p << endl;
							MPI_Recv(AT[k][k], kb*kb, MPI_DOUBLE, block_rank[k + k*p], 100000+(k+k*p), MPI_COMM_WORLD, &stat);
						}
                        // cout << "Recv A[" << k << "][" << k << "] on rank=" << myrank << endl;
					}
                } // GEQRT Recv task End

                for (int j=k+1; j<p; j++)
                {
                    const int jb = min(m-j*nb,nb); // tile col size

                    if (block_rank[k + j*p] == myrank)
                    {
                        #pragma omp task depend(in:AT[k][k],TT[k][k]) depend(inout:AT[k][j])
                        LARFB(PlasmaLeft, PlasmaTrans, 
                            kb, jb, kb, bs, 
                            AT[k][k], kb, TT[k][k], bs, AT[k][j], kb);
                        
                        cout << "LARFB(" << k << "," << j << ") on rank=" << myrank << endl;
                    }

                    if (block_rank[k + j*p] == myrank)
                    {
                        // LARFB Send task
                        // Send A[k][j] to i-direction
                    }

                    if (block_rank[k + j*p] != myrank)
                    {
                        // LARFB Recv task
                        // Recv A[k][j] from the ranks of tile A[k][j]
                    }
                }

                for (int i=k+1; i<p; i++)
                {
                    const int ib = min(m-i*nb,nb); // tile row size

                    if (block_rank[i + k*p] == myrank)
                    {
                        #pragma omp task depend(inout:AT[k][k],AT[i][k]) depend(out:TT[i][k])
                        TSQRT(kb, kb, ib, kb, bs, 
                            AT[k][k], kb, AT[i][k], ib, TT[i][k], bs);
                        
                        cout << "TSQRT(" << i << "," << k << ") on rank=" << myrank << endl;
                    }

                    if (block_rank[i + k*p] == myrank)
                    {
                        // TSQRT Send task
                        // Send A[i][k] and T[i][k] to j-direction
                    }

                    if (block_rank[i + k*p] != myrank)
                    {
                        // TSQRT Recv task
                        // Recv A[i][k] and T[i][k] from the ranks of tile A[i][k]
                    }
                    
                    for (int j=k+1; j<p; j++)
                    {
                        const int jb = min(m-j*nb,nb); // tile col size

                        if (block_rank[i + j*p] == myrank)
                        {
                            #pragma omp task depend(in:AT[i][k],TT[i][k]) depend(inout:AT[k][j],AT[i][j])
                            SSRFB(PlasmaLeft, PlasmaTrans, 
                                kb, jb, ib, jb, kb, bs,
                                AT[i][k], ib, TT[i][k], bs, AT[k][j], kb, AT[i][j], ib);
                            
                            cout << "SSRFB(" << i << "," << j << ") on rank=" << myrank << endl;
                        }

                        if (block_rank[i + j*p] == myrank)
                        {
                            // SSRFB Send task
                            // Send A[k][j] to i-direction
                        }

                        if (block_rank[i + j*p] != myrank)
                        {
                            // SSRFB Recv task
                            // Recv A[k][j] from the ranks of tile A[i][j]
                        }
                    }
                }
            }
        }
    }
    // end of main loop
    //////////////////////////////////////////////////////////////////////

    // Timer stop
    if (myrank == 0)
    {
        time = omp_get_wtime() - time;
        cout << m << ", " << nb << ", " << bs << ", " << time << endl;
    }

	for (int j=0; j<p; j++)
	{
		const int jb = min(m-j*nb,nb);

		for (int i=0; i<p; i++)
		{
			const int ib = min(m-i*nb,nb);

			if (block_rank[i + j*p] == myrank)
            {
				MPI_Bcast(AT[i][j], ib*jb, MPI_DOUBLE, myrank, MPI_COMM_WORLD);
                MPI_Bcast(TT[i][j], bs*jb, MPI_DOUBLE, myrank, MPI_COMM_WORLD);
            }
			else
            {
				MPI_Bcast(AT[i][j], ib*jb, MPI_DOUBLE, block_rank[i + j*p], MPI_COMM_WORLD);
                MPI_Bcast(TT[i][j], bs*jb, MPI_DOUBLE, block_rank[i + j*p], MPI_COMM_WORLD);
            }
		}
	}

    //////////////////////////////////////////////////////////////////////
    #ifdef DEBUG
    if (myrank == 0)
    {
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
                delete [] QT[i][j];
    }
    #endif
    //////////////////////////////////////////////////////////////////////

	MPI_Finalize();

	delete [] block_rank;

    for (int j=0; j<p; j++)
        for (int i=0; i<p; i++)
        {
            delete [] AT[i][j];
            delete [] TT[i][j];
        }
    delete [] A;

    return EXIT_SUCCESS;
}
