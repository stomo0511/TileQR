//
//  StaticPipeline
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
	// Static Pipeline tile QR
	#pragma omp parallel firstprivate(ttime)
	{
		int tk = 0;
		int tj = omp_get_thread_num();

		#ifdef TRAC
		double start_t;
		#endif

		while (tj >= NT)
		{
			tk++;
			tj = tj - NT + tk;
		}

		int ti = tk;
		
		while (tk < min(MT,NT) && tj < NT)
		{
			int next_i = ti;
			int next_j = tj;
			int next_k = tk;
			
			next_i++;
			if (next_i == MT)
			{
				next_j += omp_get_num_threads();

				while (next_j >= NT && next_k < min(MT,NT))
				{
					next_k++;
					next_j = next_j - NT + next_k;
				}
				next_i = next_k;
			}
			
			if (tj == tk)
			{
				if (ti == tk)
				{
					// GEQRT
					if (tk != 0)
						Pt.check_waitIJK( tk, tk, tk-1 );	// Check for SSRFB_(tk,tk,tk-1)

					{
						#ifdef TRAC
						start_t = omp_get_wtime();
						#endif

						GEQRT( A(tk,tk), T(tk,tk) );

						#ifdef TRAC
						#pragma omp critical
						cout << omp_get_thread_num() << ", 0, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << tk << "," << tk << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
						#endif
					}

					#ifdef DEBUG
					#pragma omp critical
					cout << "GEQRT(" << tk << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif

					Pt.setIJK(tk, tk, tk, DONE);			// Progress table update
				}   // GEQRT END
				else {
					// TSQRT
					if (tk != 0)
						Pt.check_waitIJK( ti, tk, tk-1 );	// Check for SSRFB_(ti,tk,tk-1)
					
					{
						#ifdef TRAC
						start_t = omp_get_wtime();
						#endif

						TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

						#ifdef TRAC
						#pragma omp critical
						cout << omp_get_thread_num() << ", 1, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << ti << "," << tk << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
						#endif
					}
					
					#ifdef DEBUG
					#pragma omp critical
					cout << "TSQRT(" << ti << "," << tk << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif

					Pt.setIJK(ti, tk, tk, DONE);			// Progress table update
				}   // TSQRT END
			} // Decomposition Kernel END
			else // Update Kernel
			{
				if (ti == tk)
				{
					// LARFB
					Pt.check_waitIJK( tk, tk, tk );			// Check for GEQRT_(tk,tk,tk)
					if (tk != 0)
						Pt.check_waitIJK( tk, tj, tk-1 );	// Check for SSRFB_(tk,tj,tk-1)
					
					{
						#ifdef TRAC
						start_t = omp_get_wtime();
						#endif

						LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );

						#ifdef TRAC
						#pragma omp critical
						cout << omp_get_thread_num() << ", 2, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << tk << "," << tj << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
						#endif
					}

					#ifdef DEBUG
					#pragma omp critical
					cout << "LARFB(" << tk << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif

				}   // LARFB END
				else
				{
					// SSRFB
					Pt.check_waitIJK( ti, tk, tk );			// Check for TSQRT_(ti,tk,tk)
					if (tk != 0)
						Pt.check_waitIJK( ti, tj, tk-1 );	// Check for SSRFB_(ti,tj,tk-1)

					{
						#ifdef TRAC
						start_t = omp_get_wtime();
						#endif

						SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );

						#ifdef TRAC
						#pragma omp critical
						cout << omp_get_thread_num() << ", 3, " << start_t - ttime << ", " << omp_get_wtime() - ttime << ", (" << ti << "," << tj << "," << tk << "), " << omp_get_wtime() - start_t << "\n";
						#endif
					}

					#ifdef DEBUG
					#pragma omp critical
					cout << "SSRFB(" << ti << "," << tj << "," << tk << ") : " << omp_get_thread_num() << " : " << omp_get_wtime() - ttime << "\n";
					#endif

					Pt.setIJK(ti, tj, tk, DONE);			// Progress table update
				}    // SSRFB END
			} // Update Kernel END

			ti = next_i;
			tj = next_j;
			tk = next_k;
		} // while-LOOP end
	} // End of outer most loop
	// Static Pipeline QR END
	////////////////////////////////////////////////////////////////////////////
}
