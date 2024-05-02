//=============================================================================
// Name        : BlockChol.cpp
// MPI Blocking communication, Right looking version
// Author      : Tomohiro Suzuki
// Version     : 2023/01/31
//=============================================================================

#include <iostream>
#include <iomanip>
#include <cstdlib>
#include <cassert>
#include <cmath>
#include <omp.h>
#include <mpi.h>
#include "mpi-detach.h"

#ifndef MKL
	#include <cblas.h>
	#include <lapacke.h>
#else
	#include <mkl.h>
#endif

using namespace std;

#ifdef TRACE
#include "trace.h"
#endif

// Generate random LOWER tiled matrix
void Gen_rand_lower_tiled_mat(const int i, const int j, const int m, const int n, const int mb, const int nb, double* At)
{
	assert(i>=j);

    const int p = (m % mb == 0) ? m/mb : m/mb+1;
    const int q = (n % nb == 0) ? n/nb : n/nb+1;
    const int ib = min(m-i*mb,mb);
    const int jb = min(n-j*nb,nb);

    if (i == j) // 対角タイル                                                         
        for (int jj=0; jj<jb; jj++)
            for (int ii=0; ii<ib; ii++)
                if (ii == jj) // 対角要素                                 
                    At[ii+jj*ib] = (double)(m);
                else if (ii > jj)
                    At[ii+jj*ib] = 1.0 - 2*(double)rand() / RAND_MAX;
                else
                    At[ii+jj*ib] = 0.0;
    else if (i > j) // 下三角タイル                                                   
        for (int ii=0; ii<ib*jb; ii++)
            At[ii] = 1.0 - 2*(double)rand() / RAND_MAX;
}

///////////////////////////////////////////////////////////////////////////////
// Set global variables
int myrank;  // myrank
int nprocs;  // number of MPI process
int nthred;  // number of threads

///////////////////////////////////////////////////////////////////////////////
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

int main(int argc, char **argv)
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

	// Usage "a.out [size of matrix: m ] [block width]"
	assert(argc > 2);

	const int m = atoi(argv[1]);     // # rows and columns <- the matrix is square
	const int lda = m;
	const int nb = atoi(argv[2]);    // Block size
	assert(nb <= m);

	// Number of blocks
	const int p = (m % nb == 0) ? m/nb : m/nb+1;

	// Set block rank
	int* block_rank = new int [p*p];
	get_block_rank(block_rank, p);

	///////////////////////////////////////////////////////////////////////////
	// Allocate arrays
	double* T[p][p];    // Tiled matrix
	for (int j=0; j<p; j++)
	{
		int jb = min(m-j*nb,nb);
		for (int i=j; i<p; i++)   // Lower part only
		{
			int ib = min(m-i*nb,nb);
			T[i][j] = new double [ib*jb];
			// Randomize elements of orig. matrix on each rank
			Gen_rand_lower_tiled_mat(i, j, m, m, nb, nb, T[i][j]);
		}
	}

	////////// Debug mode //////////
	#ifdef DEBUG
	// Copy T[] to OT[][] on each MPI process
	double *OT[p][p];
    #pragma omp parallel for
    for (int j=0; j<p; j++)
    {
        int jb = min(m-j*nb,nb);
        for (int i=j; i<p; i++)
        {
            int ib = min(m-i*nb,nb);
            OT[i][j] = new double [ib*jb];
			cblas_dcopy(ib*jb, T[i][j], 1, OT[i][j], 1);
        }
    }
    #endif
	
	#ifdef COMM
	int lsend = 0, lrecv = 0;
	#endif

	/////////////////////////////////////////////////////////////////////////
	// Cholesky start
	double timer = omp_get_wtime();    // Timer start

	#pragma omp parallel
	{
		#pragma omp master
		{
			for (int k=0; k<p; k++)
			{
				int kb = min(m-k*nb,nb);

				///////////////////////////////////////////////////////////////
				// POTRF
				if (block_rank[k + k*p] == myrank)
				{
					#pragma omp task depend(inout:T[k][k])
					{
						#ifdef TRACE
						trace_cpu_start();
						trace_label("Red", "dpotrf");
						#endif

						assert(0 == LAPACKE_dpotrf(LAPACK_COL_MAJOR, 'L', kb, T[k][k], kb));

						#ifdef TRACE
						trace_cpu_stop("Red");
						#endif

					}
				}

				///////////////////////////////////////////////////////////////
				// POTRF Send task
				if (block_rank[k + k*p] == myrank)
				{
					for (int dst=0; dst<nprocs; dst++)
					{
						int send_flag = 0;
						for (int kk=k+1; kk<p; kk++)
						{
							if (dst == block_rank[kk + k*p])
							{
								send_flag = 1;
								break;
							}
						}
						if (send_flag && dst != myrank)
						{
							#pragma omp task depend(in: T[k][k])
							{
								#ifdef TRACE
								trace_cpu_start();
								trace_label("Plum", "Send");
								#endif

								MPI_Send(T[k][k], kb*kb, MPI_DOUBLE, dst, k+k*p, MPI_COMM_WORLD);

								#ifdef COMM
								#pragma omp critical
								lsend++;
								#endif

								#ifdef TRACE
								trace_cpu_stop("Plum");
								#endif
							}
						}
					}
				} // POTRF Send task End

				///////////////////////////////////////////////////////////////
				// POTRF Receive task
				if (block_rank[k + k*p] != myrank)
				{
					int recv_flag = 0;
					MPI_Status st;

					for (int i=k+1; i<p; i++)
					{
						if (myrank == block_rank[i + k*p])
						{
							recv_flag = 1;
							break;
						}
					}
					if (recv_flag)
					{
						omp_event_handle_t event_handle;
						#pragma omp task detach(event_handle) depend(out: T[k][k])
						{
							MPI_Request req;

							#ifdef TRACE
							trace_cpu_start();
							trace_label("Gray", "Receive");
							#endif

							MPI_Irecv(T[k][k], kb*kb, MPI_DOUBLE, block_rank[k + k*p], k+k*p, MPI_COMM_WORLD, &req);
							MPIX_Detach(&req, (MPIX_Detach_callback *)omp_fulfill_event, (void*)event_handle);

							#ifdef COMM
							#pragma omp critical
							lrecv++;
							#endif

							#ifdef TRACE
							trace_cpu_stop("Gray");
							#endif
						}
						recv_flag = 0;
					}
				} // POTRF Receive task End

				for (int i=k+1; i<p; i++)
				{
					int ib = min(m-i*nb,nb);

					///////////////////////////////////////////////////////
					// TRSM
					if (block_rank[i + k*p] == myrank)
					{
						#pragma omp task depend(in:T[k][k]) depend(inout:T[i][k])
						{
							#ifdef TRACE
							trace_cpu_start();
							trace_label("Green", "dtrsm");
							#endif

							cblas_dtrsm(CblasColMajor, CblasRight, CblasLower, CblasTrans, CblasNonUnit, 
								ib, kb, 1.0, T[k][k], kb, T[i][k], ib);

							#ifdef TRACE
							trace_cpu_stop("Green");
							#endif
						}
					}

					///////////////////////////////////////////////////////////////
					// TRSM Send task
					if (block_rank[i + k*p] == myrank)
					{
						int* send_flags = new int[nprocs];
						for (int ii=0; ii<nprocs; ii++)
							send_flags[ii] = 0;
							
						for (int ii = k + 1; ii < i; ii++)
						{
							if (!send_flags[block_rank[i + ii*p]])
								send_flags[block_rank[i + ii*p]] = 1;
						}
						for (int ii = i + 1; ii < p; ii++)
						{
							if (!send_flags[block_rank[ii + i*p]])
								send_flags[block_rank[ii + i*p]] = 1;
						}
						if (!send_flags[block_rank[i + i*p]])
							send_flags[block_rank[i + i*p]] = 1;

						for (int dst=0; dst<nprocs; dst++)
						{
							if (send_flags[dst] && dst != myrank)
							{
								#pragma omp task depend(in: T[i][k])
								{
									#ifdef TRACE
									trace_cpu_start();
									trace_label("Plum", "Send");
									#endif

									MPI_Send(T[i][k], ib*kb, MPI_DOUBLE, dst, i+k*p, MPI_COMM_WORLD);

									#ifdef COMM
									#pragma omp critical
									lsend++;
									#endif

									#ifdef TRACE
									trace_cpu_stop("Plum");
									#endif
								}
							}
						}
						delete [] send_flags;
					} // TRSM Send task End

					///////////////////////////////////////////////////////////////
					// TRSM Receive task
					if (block_rank[i + k*p] != myrank)
					{
						int recv_flag = 0;
						MPI_Status st;

						for (int ii = k + 1; ii < i; ii++)
						{
							if (block_rank[i + ii*p] == myrank)
								recv_flag = 1;
						}
						for (int ii = i + 1; ii < p; ii++)
						{
							if (block_rank[ii + i*p] == myrank)
								recv_flag = 1;
						}
						if (block_rank[i + i*p] == myrank)
							recv_flag = 1;

						if (recv_flag)
						{
							omp_event_handle_t event_handle;
							#pragma omp task detach(event_handle) depend(out: T[i][k])
							{
								MPI_Request req;

								#ifdef TRACE
								trace_cpu_start();
								trace_label("Gray", "Receive");
								#endif

								MPI_Irecv(T[i][k], ib*kb, MPI_DOUBLE, block_rank[i + k*p], i+k*p, MPI_COMM_WORLD, &req);
								MPIX_Detach(&req, (MPIX_Detach_callback *)omp_fulfill_event, (void*)event_handle);

								#ifdef COMM
								#pragma omp critical
								lrecv++;
								#endif

								#ifdef TRACE
								trace_cpu_stop("Gray");
								#endif
							}
							recv_flag = 0;
						}
					} // TRSM Receive task End

					///////////////////////////////////////////////////////////////
					// SYRK
					if (block_rank[i + i*p] == myrank)
					{
						#pragma omp task depend(in:T[i][k]) depend(inout:T[i][i])
						{
							#ifdef TRACE
							trace_cpu_start();
							trace_label("Cyan", "dsyrk");
							#endif

							cblas_dsyrk(CblasColMajor, CblasLower, CblasNoTrans, 
								ib, kb, -1.0, T[i][k], ib, 1.0, T[i][i], ib);

							#ifdef TRACE
							trace_cpu_stop("Cyan");
							#endif
						}
					}

					for (int j=k+1; j<i; j++)
					{
						int jb = min(m-j*nb,nb);

						///////////////////////////////////////////////////
						// GEMM
						if (block_rank[i + j*p] == myrank)
						{
							#pragma omp task depend(in:T[i][k], T[j][k]) depend(inout:T[i][j])
							{
								#ifdef TRACE
								trace_cpu_start();
								trace_label("Blue", "dgemm");
								#endif

								cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans, ib, jb, kb, 
									-1.0, T[i][k], ib, T[j][k], nb, 1.0, T[i][j], ib);

								#ifdef TRACE
								trace_cpu_stop("Blue");
								#endif
							}
						}
					} // End of J-loop
				} // End of I-loop
			} // End of K-loop
			#pragma omp taskwait
		} // End of single block
	} // End of parallel block

	timer = omp_get_wtime() - timer;   // Timer stop
	// Cholesky end

	for (int j=0; j<p; j++)
	{
		int jb = min(m-j*nb,nb);
		for (int i=j; i<p; i++)
		{
			int ib = min(m-i*nb,nb);
			if (block_rank[i + j*p] == myrank)
				MPI_Bcast(T[i][j], ib*jb, MPI_DOUBLE, myrank, MPI_COMM_WORLD);
			else
				MPI_Bcast(T[i][j], ib*jb, MPI_DOUBLE, block_rank[i + j*p], MPI_COMM_WORLD);
		}
	}

	if (myrank == 0)
	{
		cout << nprocs << ", ";
		cout << omp_get_max_threads() << ", ";
		cout << m << ", " << nb << ", " << timer;

	}

	#ifdef COMM
	int gsend = 0, grecv = 0;
	MPI_Reduce( &lsend, &gsend, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );
	MPI_Reduce( &lrecv, &grecv, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD );

	if (myrank == 0)
	{
		// assert(gsend == grecv);
		cout << ", " << gsend << ", " << grecv;
	}
	#endif

	if (myrank == 0)
		cout << endl;

	////////// Debug mode //////////
	#ifdef DEBUG
	double norm = 0.0;

	for (int j=0; j<p; j++)
	{
		int jb = min(m-j*nb,nb);
		for (int i=j; i<p; i++)
		{
			int ib = min(m-i*nb,nb);
			if (block_rank[i + j*p] == myrank)
			{
				for (int k=0; k<=j; k++)
				{
					int kb = min(m-k*nb,nb);

					double alpha = -1.0;
					double beta = 1.0;
					cblas_dgemm(CblasColMajor, CblasNoTrans, CblasTrans,
						ib, jb, kb, alpha, T[i][k], ib, T[j][k], jb,
						beta, OT[i][j], ib);

					if ((j == k) && (i == k))
					{
						for (int jj=0; jj<jb; jj++)
							for (int ii=0; ii<jj; ii++)
								(OT[i][j])[ii + jj*ib] = 0.0;
					}
				}
				norm += cblas_dnrm2(ib*jb, OT[i][j], 1);
			}
		}
	}

	double gnorm = 0.0;
	MPI_Reduce( &norm, &gnorm, 1, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD );

	if (myrank == 0)
		cout << "|| A - LL^T ||_2 = " << gnorm << endl;

	for (int j=0; j<p; j++)
		for (int i=j; i<p; i++)
			delete [] OT[i][j];
	#endif
	////////// Debug mode //////////

	MPI_Finalize();

	///////////////////////////////////////////////////////////////////////////
	// Free arrays
	delete [] block_rank;

	for (int j=0; j<p; j++)
		for (int i=j; i<p; i++)
			delete [] T[i][j];

	return EXIT_SUCCESS;
}
