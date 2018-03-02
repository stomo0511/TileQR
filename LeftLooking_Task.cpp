//
//  LeftLooking Dynamic Scheduring version
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <omp.h>

#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>

//#define DEBUG

using namespace std;

void tileQR( const int MT, const int NT, TMatrix& A, TMatrix& T )
{
	// Progress table
	int **Ap;

	Ap = (int **)malloc( sizeof(int*) * MT);
	for (int i=0; i<MT; i++)
		Ap[i] = (int *)malloc( sizeof(int) * NT);

	double ttime = omp_get_wtime();

	////////////////////////////////////////////////////////////////////////////
	// Left Looking tile QR
	#pragma omp parallel firstprivate(ttime)
	{
		#pragma omp single
		{
			for (int tk = 0; tk < NT; tk++)
			{
				for (int tl = 0; tl < min(MT,tk); tl++)
				{
					// LARFB
					#pragma omp task depend(in:Ap[tl][tl]) depend(inout:Ap[tl][tk])
					{
						LARFB( PlasmaLeft, PlasmaTrans, A(tl,tl), T(tl,tl), A(tl,tk) );

						#ifdef DEBUG
						#pragma omp critical
						cout << "LARFB(" << tl << "," << tk << "," << tl << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
						#endif
					}

					for (int ti = tl+1; ti < MT; ti++)
					{
						// SSRFB
						#pragma omp task depend(in:Ap[ti][tl]) depend(inout:Ap[tl][tk], Ap[ti][tk])
						{
							SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tl), T(ti,tl), A(tl,tk), A(ti,tk) );

							#ifdef DEBUG
							#pragma omp critical
							cout << "SSRFB(" << ti << "," << tk << "," << tl << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
							#endif
						}
					} // End of i-loop
				} // End of l-loop

				if ( tk < MT )
				{
					// GEQRT
					#pragma omp task depend(inout:Ap[tk][tk])
					{
						GEQRT( A(tk,tk), T(tk,tk) );

						#ifdef DEBUG
						#pragma omp critical
						cout << "GEQRT(" << tk << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
						#endif
					}

					for (int ti = tk+1; ti < MT; ti++)
					{
						// TSQRT
						#pragma omp task depend(inout:Ap[tk][tk]) depend(out:Ap[ti][tk])
						{
							TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

							#ifdef DEBUG
							#pragma omp critical
							cout << "TSQRT(" << ti << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
							#endif
						}
					} // End of i-loop
				} // End if ( k < MT )
			} // End of k-loop
		} // End of single section
	} // End of parallel section
  // Left Looking tile QR END
  ////////////////////////////////////////////////////////////////////////////
}
