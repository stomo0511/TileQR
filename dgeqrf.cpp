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

#include "Check_Accuracy.h"

using namespace std;

// Generate random LOWER matrix
void Gen_rand_mat(const int m, const int n, double *A)
{
	srand(20240417);

    for (int j=0; j<n; j++)
	    for (int i=0; i<m; i++)
            A[i+j*m] = 1.0 - 2*(double)rand() / RAND_MAX;
}

// Show matrix
void Show_mat(const int m, const int n, double *A)
{
	for (int i=0; i<m; i++)
    {
		for (int j=0; j<n; j++)
			cout << A[i + j*m] << ", ";
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

    //Usage "a.out [size of matrix: M ]"
    if (argc < 1)
    {
        cout << "Usage: a.out [size of matrix: M ]" << endl;
        return EXIT_FAILURE;
    }

    const int m =  atoi(argv[1]);  // matrix size

    double *A = new double [m*m];    // Original matrix: A
    Gen_rand_mat(m,m,A);             // Randomize elements of A

    double *T = new double [m];      // Contain scalers of elementary reflectors
	
	////////// Debug mode //////////
	#ifdef DEBUG
	double *OA = new double [m*m];
	cblas_dcopy(m*m, A, 1, OA, 1);
    double *R = new double [m*m];
	#endif
	////////// Debug mode //////////

//    for (int i=0; i<m; i++)
//     {
//         for (int j=0; j<m; j++)
//                 cout << A[i+j*m] << ", ";
//         cout << endl;
//     }
//     cout << endl;

    // Timer start
    double time = omp_get_wtime();
	
    //////////////////////////////////////////////////////////////////////
    LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, m, A, m, T);
    //////////////////////////////////////////////////////////////////////
	
    // Timer stop
    time = omp_get_wtime() - time;
    cout << m << ", " << time << endl;
	
    // for (int i=0; i<m; i++)
    // {
    //     for (int j=0; j<m; j++)
    //             cout << A[i+j*m] << ", ";
    //     cout << endl;
    // }

    #ifdef DEBUG
    // Copy upper triangular part of A to R
    for (int i=0; i<m; i++)
        for (int j=0; j<m; j++)
            if (i <= j)
                R[i+j*m] = A[i+j*m];
            else
                R[i+j*m] = 0.0;

    // Regenerate orthogonal matrix Q
    LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, m, m, A, m, T);

    Check_Accuracy(m,m,OA,A,R);
    // // Check Accuracy END
    // //////////////////////////////////////////////////////////////////////

	delete [] OA;
    delete [] R;
    #endif

    delete [] A;
    delete [] T;

    return EXIT_SUCCESS;
}
