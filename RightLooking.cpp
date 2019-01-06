//
//  RightLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

//#define COUT
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
			#pragma omp single
			{
				GEQRT( A(tk,tk), T(tk,tk) );

				#ifdef COUT
				cout << "GEQRT(" << tk << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
				#endif
				#ifdef ANIM
				cout << "GL," << tk << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
				#endif
			}

#pragma omp for
			for (int tj=tk+1; tj < NT; tj++)
			{
				LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );

				#ifdef COUT
				#pragma omp critical
				cout << "LARFB(" << tk << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
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
					TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

					#ifdef COUT
					cout << "TSQRT(" << ti << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif
					#ifdef ANIM
					#pragma omp critical
					cout << "TL," << ti << "," << tk << "," << tk << "," << omp_get_wtime() - ttime << endl;
					#endif
				}

#pragma omp for
				for (int tj=tk+1; tj < NT; tj++)
				{
					SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );

					#ifdef COUT
					#pragma omp critical
					cout << "SSRFB(" << ti << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
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
