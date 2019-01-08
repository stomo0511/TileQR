//
//  LeftLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#define TRAC

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <omp.h>

#include "Progress.hpp"
#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>

void tileQR( const int MT, const int NT, TMatrix& A, TMatrix& T )
{
	// Progress Table
	Progress_Table Pt( MT, NT, min(MT,NT) );
	Pt.Init();
	Pt.setIJK(0, 0, 0, NYET);
	
	double ttime = omp_get_wtime();

	////////////////////////////////////////////////////////////////////////////
	// Left Looking tile QR
	#pragma omp parallel for schedule(static,1) firstprivate(ttime)
	for (int tk = 0; tk < NT; tk++)
	{
		#ifdef TRAC
		double start_t;
		#endif
			
		for (int tl = 0; tl < min(MT,tk); tl++)
		{

			Pt.check_waitIJK(tl, tl, tl);	// Check for GEQRT_(tl,tl,tl)

			// LARFB
			{
				#ifdef TRAC
				start_t = omp_get_wtime();
				#endif

				LARFB( PlasmaLeft, PlasmaTrans, A(tl,tl), T(tl,tl), A(tl,tk) );

				#ifdef DEBUG
				#pragma omp critical
				cout << "LARFB(" << tl << "," << tk << "," << tl << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
				#endif
				#ifdef TRAC
				#pragma omp critical
				cout << omp_get_thread_num() << ", 2, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << tl << "," << tk << "," << tl << "), " << omp_get_wtime() - start_t << "\n";
				#endif
			}

			for (int ti = tl+1; ti < MT; ti++)
			{
				Pt.check_waitIJK(ti, tl, tl);	// Check for TSQRT_(ti,tl,tl)

				// SSRFB
				{
					#ifdef TRAC
					start_t = omp_get_wtime();
					#endif

					SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tl), T(ti,tl), A(tl,tk), A(ti,tk) );

					#ifdef DEBUG
					#pragma omp critical
					cout << "SSRFB(" << ti << "," << tk << "," << tl << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif
					#ifdef TRAC
					#pragma omp critical
					cout << omp_get_thread_num() << ", 3, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << ti << "," << tk << "," << tl << "), " << omp_get_wtime() - start_t << "\n";
					#endif
				}
			} // End of i-loop
		} // End of l-loop

		if ( tk < MT )
		{
			// GEQRT
			{
				#ifdef TRAC
				start_t = omp_get_wtime();
				#endif

				GEQRT( A(tk,tk), T(tk,tk) );

				#ifdef DEBUG
				#pragma omp critical
				cout << "GEQRT(" << tk << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
				#endif
				#ifdef TRAC
				#pragma omp critical
				cout << omp_get_thread_num() << ", 0, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << tk << "," << tk << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
				#endif
			}

			// Progress table update
			Pt.setIJK(tk, tk, tk, DONE);
				
			for (int ti = tk+1; ti < MT; ti++)
			{
				// TSQRT
				{
					#ifdef TRAC
					start_t = omp_get_wtime();
					#endif

					TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

					#ifdef DEBUG
					#pragma omp critical
					cout << "TSQRT(" << ti << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif
					#ifdef TRAC
					#pragma omp critical
					cout << omp_get_thread_num() << ", 1, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << ti << "," << tk << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
					#endif
				}

				// Progress table update
				Pt.setIJK(ti, tk, tk, DONE);
					
			} // End of i-loop
		} // End if ( k < MT )
	} // End of k-loop
  // Left Looking tile QR END
  ////////////////////////////////////////////////////////////////////////////
}
