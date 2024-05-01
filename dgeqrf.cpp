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

using namespace std;

// Generate random LOWER matrix
void Gen_rand_mat(const long int seed, const int m, const int n, double *A)
{
	srand(seed);

    for (int j=0; j<n; j++)
	    for (int i=0; i<m; i++)
            A[i+j*m] = 1.0 - 2*(double)rand() / RAND_MAX;
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

void Check_Accuracy( const int m, const int n, double *mA, double *mQ, double *mR )
{
    //////////////////////////////////////////////////////////////////////
    // Check Orthogonarity

    // Set Id to the identity matrix
    int mn = min(m,n);

    // Id: identity matrix
    double* Id = new double[ mn * mn ];
    LAPACKE_dlaset(LAPACK_COL_MAJOR, 'g', mn, mn, 0.0, 1.0, Id, mn);

    double alpha = 1.0;
    double beta  = -1.0;
        
    cblas_dsyrk(CblasColMajor, CblasUpper, CblasTrans, 
            n, m, alpha, mQ, m, beta, Id, n);
        
    double ortho = LAPACKE_dlansy(LAPACK_COL_MAJOR, 'I', 'U', 
                        mn, Id, mn);

    // normalize the result
    // |Id - Q^T * Q|_oo / n
    // ortho /= mn;
    std::cout << "norm(I-Q'*Q) = " << ortho << std::endl;

    // Check Orthogonarity END
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Check Residure

    // |A|_oo
    double normA = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', m, n, mA, m);

    alpha = -1.0;
    beta  =  1.0;
    cblas_dgemm(CblasColMajor, CblasNoTrans, CblasNoTrans, 
            m, n, m, alpha, mQ, m, mR, m, beta, mA, m);
        
    double normQ = LAPACKE_dlange(LAPACK_COL_MAJOR, 'I', 
                    m, n, mA, m);

    // normalize the result
    // |A-QR|_oo / (|A|_oo * n)
    // normQ /= (normA * N);
    std::cout << "norm(A-Q*R) = " << normQ << std::endl;

    // Check Residure END
    //////////////////////////////////////////////////////////////////////

    delete [] Id;
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

    const int m =  atoi(argv[1]);    // matrix size

    double *A = new double [m*m];    // Original matrix: A
    Gen_rand_mat(20240417, m, m, A); // Randomize elements of A

    double *T = new double [m];      // Contain scalers of elementary reflectors
	
	////////// Debug mode //////////
	#ifdef DEBUG
	double *OA = new double [m*m];
	cblas_dcopy(m*m, A, 1, OA, 1);
    double *R = new double [m*m];
	#endif
	////////// Debug mode //////////

    // Timer start
    double time = omp_get_wtime();
	
    LAPACKE_dgeqrf(LAPACK_COL_MAJOR, m, m, A, m, T);
	
    // Timer stop
    time = omp_get_wtime() - time;
    cout << m << ", " << time << endl;
	
    #ifdef DEBUG
    // Copy upper triangular part of A to R
    for (int i=0; i<m; i++)
        for (int j=0; j<m; j++)
            if (i <= j)
                R[i+j*m] = A[i+j*m];
            else
                R[i+j*m] = 0.0;

    // Regenerate orthogonal matrix Q
    assert(EXIT_SUCCESS == LAPACKE_dorgqr(LAPACK_COL_MAJOR, m, m, m, A, m, T));

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
