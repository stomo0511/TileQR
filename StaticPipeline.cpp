//
//  StaticPipeline
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <omp.h>

#include "Progress.hpp"
#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>

#ifndef __Test__Min__
#define __Test__Min__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __Test__Min__

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

					GEQRT( A(tk,tk), T(tk,tk) );

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
					
					TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );
					
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
					
					LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );

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

					SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );

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
