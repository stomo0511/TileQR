//
//  RightLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

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
	#ifdef ANIM
	cout << "Kernel,Ii,Ij,Ik,Time\n";
	#endif

	double ttime = omp_get_wtime();

	for (int tk=0; tk < min(MT,NT); tk++ )
	{
		#pragma omp parallel firstprivate(ttime)
		{
			#ifdef TRAC
			double start_t;
			#endif

			#pragma omp single
			{
				#ifdef TRAC
				start_t = omp_get_wtime();
				#endif

				GEQRT( A(tk,tk), T(tk,tk) );

				#ifdef COUT
				cout << "GEQRT(" << tk << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
				#endif
				#ifdef TRAC
				#pragma omp critical
				cout << omp_get_thread_num() << ", 0, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << tk << "," << tk << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
				#endif
				#ifdef ANIM
				cout << "GL," << tk << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
				#endif
			}

#pragma omp for
			for (int tj=tk+1; tj < NT; tj++)
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
				cout << omp_get_thread_num() << ", 2, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << tk << "," << tj << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
				#endif
				#ifdef ANIM
				#pragma omp critical
				cout << "LL," << tk << "," << tj << "," << tk << "," << omp_get_wtime() - ttime << endl;
				#endif
			} // j-LOOP END

			for (int ti=tk+1; ti < MT; ti++)
			{
				#pragma omp single
				{
					#ifdef TRAC
					start_t = omp_get_wtime();
					#endif

					TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

					#ifdef COUT
					cout << "TSQRT(" << ti << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif
					#ifdef TRAC
					#pragma omp critical
					cout << omp_get_thread_num() << ", 1, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << ti << "," << tk << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
					#endif
					#ifdef ANIM
					#pragma omp critical
					cout << "TL," << ti << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				}

#pragma omp for
				for (int tj=tk+1; tj < NT; tj++)
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
					cout << omp_get_thread_num() << ", 3, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << ti << "," << tj << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
					#endif
					#ifdef ANIM
					#pragma omp critical
					cout << "SL," << ti << "," << tj << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				} // j-LOOP END
			} // i-LOOP END
		} // parallel section END
	} // k-LOOP END
	// Right Looking tile QR END
	//////////////////////////////////////////////////////////////////////
}
