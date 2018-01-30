/*
 * geqrf.cpp
 *
 *  Created on: 2016/07/21
 *      Author: stomo
 */


#include <iostream>
#include <omp.h>
#include <cassert>
#include <cstdlib>

#include "Check_Accuracy.hpp"

#include <mkl_lapacke.h>

using namespace std;

int main(int argc, const char * argv[])
{
	if (argc < 3)
	{
		cerr << "Usage: a.out [M] [N]\n";
		exit (1);
	}

	const int M =  atoi(argv[1]);  // n. of rows of the matrix
	const int N =  atoi(argv[2]);  // n. of columns of the matrix

	assert( M >= N );

	#ifdef DEBUG
	cout << "M = " << M << ", N = " << N << endl;
	cout << "clock precision = " << omp_get_wtick() << endl;
	#endif

	lapack_int info;
	unsigned seed = 20140105;

	double *A = new double [ M*N ];
	for (int i=0; i<M*N; i++)
		A[i] = (double)rand() / RAND_MAX;

	double *tau = new double [ M ];

	// Timer start
	double time = omp_get_wtime();

	info = LAPACKE_dgeqrf(LAPACK_COL_MAJOR, M, N, A, M, tau);
	if (info != 0)
	{
		cerr << "dgeqrf failed. info = " << info << "\n";
		return EXIT_FAILURE;
	}

	// Timer stop
	time = omp_get_wtime() - time;
	cout << M << ", " << N << ", " << time << endl;

	delete [] A;
	delete [] tau;

	return EXIT_SUCCESS;
}
