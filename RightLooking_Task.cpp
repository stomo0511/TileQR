/*
 * RightLooking_Task.cpp
 *
 *  Created on: 2016/03/31
 *      Author: stomo
 */

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <omp.h>

#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>

#ifndef __Test__Min__
#define __Test__Min__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __Test__Min__

using namespace std;

void tileQR( const int MT, const int NT, TMatrix& A, TMatrix& T )
{
	//////////////////////////////////////////////////////////////////////
	// List Item
	int **Ap;

	Ap = (int **)malloc(sizeof(int *) * MT);
	for (int i=0; i<MT; i++)
		Ap[i] = (int *)malloc(sizeof(int) * NT);

	double ttime = omp_get_wtime();

	//////////////////////////////////////////////////////////////////////
	// Right Looking tile QR Task version
	#pragma omp parallel firstprivate(ttime)
	{
		#pragma omp single
		{
			for (int tk=0; tk < min(MT,NT); tk++ )
			{
				{
					#pragma omp task depend(inout:Ap[tk][tk])
					{
						GEQRT( A(tk,tk), T(tk,tk) );

						#ifdef DEBUG
						#pragma omp critical
						cout << "GEQRT(" << tk << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
						#endif
					}
				}

				for (int tj=tk+1; tj < NT; tj++)
				{
					#pragma omp task depend(in:Ap[tk][tk]) depend(inout:Ap[tk][tj])
					{
						LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );

						#ifdef DEBUG
						#pragma omp critical
						cout << "LARFB(" << tk << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
						#endif
					}
				}

				for (int ti=tk+1; ti < MT; ti++)
				{
					{
						#pragma omp task depend(inout:Ap[tk][tk]) depend(out:Ap[ti][tk])
						{
							TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

							#ifdef DEBUG
							#pragma omp critical
							cout << "TSQRT(" << ti << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
							#endif
						}
					}

					for (int tj=tk+1; tj < NT; tj++)
					{
						#pragma omp task depend(in:Ap[ti][tk]) depend(inout:Ap[tk][tj],Ap[ti][tj])
						{
							SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );

							#ifdef DEBUG
							#pragma omp critical
							cout << "SSRFB(" << ti << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
							#endif
						}
					} // j-LOOP END
				} // i-LOOP END
			} // k-LOOP END
		} // parallel section END
	}
	// Right Looking tile QR END
	//////////////////////////////////////////////////////////////////////

	for (int i=0; i<MT; i++)
		free(Ap[i]);
	free(Ap);
}
