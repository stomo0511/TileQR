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

#include "Check_Accuracy.h"
#include "CoreBlas.h"

using namespace std;

// Generate random tiled matrix
void Gen_rand_tiled_mat(const int i, const int j, const int m, const int n, const int mb, const int nb, double *At)
{
	srand(20240417);

	const int ib = min(m-i*mb,mb);     // tile row size
	const int jb = min(n-j*nb,nb);     // tile col size

    for (int jj=0; jj<jb; jj++)
        for (int ii=0; ii<ib; ii++)
            At[ii+jj*ib] = 1.0 - 2*(double)rand() / RAND_MAX;
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
    const int n = m;
    const int mb = atoi(argv[2]);   // tile size
    const int nb = mb;
    const int ib = atoi(argv[3]);   // inner blocking size

   	const int nmb = (m % mb == 0) ? m/mb : m/mb+1;  // # of tile rows
	const int nnb = (n % nb == 0) ? n/nb : n/nb+1;  // # of tile cols

    assert(m >= mb);
    assert(mb >= ib);
	
    #ifdef DEBUG
    cout << "m = " << m << ", mb = " << mb << ", ib = " << ib << endl;
    cout << "# of threads = " << omp_get_max_threads() << endl;
    #endif
	
    double *AT[nmb][nnb];          // AT: Tiled matrix
    double *TT[nmb][nnb];          // TT: Tiled Householder matrix
    for (int j=0; j<nnb; j++)
    {
        const int nbj = min(n-j*nb,nb);     // tile col size
        for (int i=0; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size
            AT[i][j] = new double [mbi * nbj];
            Gen_rand_tiled_mat(i,j,m,n,mb,nb,AT[i][j]);

            TT[i][j] = new double [ib * nbj];
        }
    }

   for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
            cout << AT[0][0][i+j*m] << ", ";
        cout << endl;
    }
    cout << endl;

    // Timer start
    double time = omp_get_wtime();

    for (int k=0; k<min(nmb,nnb); k++)
    {
        const int mbk = min(m-k*mb,mb);     // tile row size

        GEQRT(mbk, mbk, ib, 
            AT[k][k], mbk, TT[k][k], ib);

        for (int j=k+1; j<nnb; j++)
        {
            const int nbj = min(n-j*nb,nb); // tile col size

            LARFB(PlasmaLeft, PlasmaTrans, 
                mbk, nbj, mbk, ib, 
                AT[k][k], mbk, TT[k][k], ib, AT[k][j], mbk);
        }

        for (int i=k+1; i<nmb; i++)
        {
            const int mbi = min(m-i*mb,mb); // tile row size

            TSQRT(mbk, mbk, mbi, mbk, ib, 
                AT[k][k], mbk, AT[i][k], mbi, TT[i][k], ib);
            
            for (int j=k+1; j<nnb; j++)
            {
                const int nbj = min(n-j*nb,nb); // tile col size

                SSRFB(PlasmaLeft, PlasmaTrans, 
                    mbk, nbj, mbi, nbj, mbk, ib,
                    AT[i][k], mbi, TT[i][k], ib, AT[k][j], mbk, AT[i][j], mbi);
            }
        }
    }

    // Timer stop
    time = omp_get_wtime() - time;
    cout << m << ", " << nb << ", " << ib << ", " << time << endl;

   for (int i=0; i<m; i++)
    {
        for (int j=0; j<n; j++)
            cout << AT[0][0][i+j*m] << ", ";
        cout << endl;
    }

    for (int j=0; j<nnb; j++)
        for (int i=0; i<nmb; i++)
        {
            delete [] AT[i][j];
            delete [] TT[i][j];
        }

    return EXIT_SUCCESS;
}
