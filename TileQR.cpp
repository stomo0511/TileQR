//
//  main.cpp
//
//  Created by T. Suzuki on 2022/12/22.
//

#include <iostream>
#include <omp.h>
#include <cassert>

#include "TMatrix.h"
#include "CoreBlas.h"
#include "Check_Accuracy.h"

using namespace std;

void tileQR( const int MT, const int NT, TMatrix< Tile<double> >& A, TMatrix< Tile<double> >& T );

int main(int argc, const char ** argv)
{

    //Usage "a.out [# of rows of matrix: m ] [# of columns of matrix: n ] [tile size] [inner block size]"
    assert(argc > 4);

    const int M =  atoi(argv[1]);  // n. of rows of the matrix
    const int N =  atoi(argv[2]);  // n. of columns of the matrix
    const int NB = atoi(argv[3]);  // tile size
    const int IB = atoi(argv[4]);  // inner blocking size

    assert( M >= N );
    assert( NB >= IB );
	
    #ifdef DEBUG
    cout << "M = " << M << ", N = " << N << ", NB = " << NB << ", IB = " << IB << endl;
    cout << "# of threads =" << omp_get_max_threads() << endl;
    #endif
	
    //////////////////////////////////////////////////////////////////////
    // Definitions and Initialize
    TMatrix< Tile<double> > A(M,N,NB,NB,IB);
	
    const int MT = A.mt();
    const int NT = A.nt();
	
    // refered in workspace.c of PLASMA
    TMatrix< Tile<double> > T(MT*IB,NT*NB,IB,NB,IB);
	
    // Initialize matrix A
    A.Set_Rnd( 20221222 );

    #ifdef DEBUG
    // Copy the elements of TMatrix class A to double array mA
    double *mA = new double [ M*N ];
    A.Array_Copy(mA);
    #endif
    // Definitions and Initialize　END
    //////////////////////////////////////////////////////////////////////

    // Timer start
    double time = omp_get_wtime();
	
    //////////////////////////////////////////////////////////////////////
    // tile QR variants
    tileQR(MT,NT,A,T);
    //////////////////////////////////////////////////////////////////////
	
    // Timer stop
    time = omp_get_wtime() - time;
    cout << M << ", " << NB << ", " << IB << ", " << time << endl;
	
    #ifdef DEBUG
    //////////////////////////////////////////////////////////////////////
    // Regenerate Q
    TMatrix< Tile<double> > Q(M,M,NB,NB,IB);
	
    // Set to the identity matrix
    Q.Set_Iden();
	
    // Make Orthogonal matrix Q
    dorgqr( A, T, Q );
    //////////////////////////////////////////////////////////////////////

    //////////////////////////////////////////////////////////////////////
    // Check Accuracy

    double *mQ = new double [ M*M ];
    double *mR = new double [ M*N ];

    Q.Array_Copy(mQ);
    A.Array_Copy(mR);

    for (int i=0; i<M; i++)
        for (int j=0; j<N; j++)
            if (i > j)
          	  mR[ i + M*j ] = 0.0;

    Check_Accuracy(M,N,mA,mQ,mR );
    // Check Accuracy END
    //////////////////////////////////////////////////////////////////////

    delete [] mA;
    delete [] mQ;
    delete [] mR;

    cout << "Done\n";
    #endif
	
    return EXIT_SUCCESS;
}
