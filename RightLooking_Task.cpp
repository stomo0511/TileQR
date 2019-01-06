/*
 * RightLooking_Task.cpp
 *
 *  Created on: 2016/03/31
 *      Author: stomo
 */

//#define COUT
#define TRAC
//#define ANIM

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <omp.h>

#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>

using namespace std;

void tileQR( const int MT, const int NT, TMatrix& A, TMatrix& T )
{
	// Progress table
	int **Ap;

	Ap = (int **)malloc( sizeof(int*) * MT);
	for (int i=0; i<MT; i++)
		Ap[i] = (int *)malloc( sizeof(int) * NT);

	#ifdef ANIM
	cout << "Kernel,Ii,Ij,Ik,Time\n";
	#endif

	double ttime = omp_get_wtime();

	//////////////////////////////////////////////////////////////////////
	// Right Looking tile QR Task version
	#pragma omp parallel firstprivate(ttime)
	{
		#ifdef TRAC
		double start_t;
		#endif

		#pragma omp single
		{
			for (int tk=0; tk < min(MT,NT); tk++ )
			{
				#pragma omp task depend(inout:Ap[tk][tk])
				{
					#ifdef TRAC
					start_t = omp_get_wtime();
					#endif

					GEQRT( A(tk,tk), T(tk,tk) );

					#ifdef COUT
					#pragma omp critical
					cout << "GEQRT(" << tk << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif
					#ifdef TRAC
					#pragma omp critical
					cout << omp_get_thread_num() << ", 0, " << start_t - ttime << ", " << omp_get_wtime() - ttime << " (" << tk << "," << tk << "," << tk << ")\n";
					#endif
					#ifdef ANIM
					cout << "GL," << tk << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				}

				for (int tj=tk+1; tj < NT; tj++)
				{
					#pragma omp task depend(in:Ap[tk][tk]) \
									 depend(inout:Ap[tk][tj])
					{
						#ifdef TRAC
						start_t = omp_get_wtime();
						#endif

						LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );

						#ifdef COUT
						#pragma omp critical
						cout << "LARFB(" << tk << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
						#endif
						#ifdef TRAC
						#pragma omp critical
						cout << omp_get_thread_num() << ", 2, " << start_t - ttime << ", " << omp_get_wtime() - ttime << " (" << tk << "," << tj << "," << tk << ")\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "LL," << tk << "," << tj << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}
				}

				for (int ti=tk+1; ti < MT; ti++)
				{
					#pragma omp task depend(inout:Ap[tk][tk]) \
									 depend(out:Ap[ti][tk])
					{
						#ifdef TRAC
						start_t = omp_get_wtime();
						#endif

						TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

						#ifdef COUT
						#pragma omp critical
						cout << "TSQRT(" << ti << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
						#endif
						#ifdef TRAC
						#pragma omp critical
						cout << omp_get_thread_num() << ", 1, " << start_t - ttime << ", " << omp_get_wtime() - ttime << " (" << ti << "," << tk << "," << tk << ")\n";
						#endif
						#ifdef ANIM
						#pragma omp critical
						cout << "TL," << ti << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
						#endif
					}

					for (int tj=tk+1; tj < NT; tj++)
					{
						#pragma omp task depend(in:Ap[ti][tk]) \
										 depend(inout:Ap[tk][tj], Ap[ti][tj])
						{
							#ifdef TRAC
							start_t = omp_get_wtime();
							#endif

							SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );

							#ifdef COUT
							#pragma omp critical
							cout << "SSRFB(" << ti << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
							#endif
							#ifdef TRAC
							#pragma omp critical
							cout << omp_get_thread_num() << ", 3, " << start_t - ttime << ", " << omp_get_wtime() - ttime << " (" << ti << "," << tj << "," << tk << ")\n";
							#endif
							#ifdef ANIM
							#pragma omp critical
							cout << "SL," << ti << "," << tj << "," << tk << "," << omp_get_wtime() - ttime << endl;
							#endif
						}
					} // j-LOOP END
				} // i-LOOP END
			} // k-LOOP END
		} // parallel section END
	}
	// Right Looking tile QR task END
	//////////////////////////////////////////////////////////////////////
}
