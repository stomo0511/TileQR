//
//  RightLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

//#define COUT
#define DEBUG

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <algorithm>
#include <omp.h>

#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>

using namespace std;

void tileTRD( const int MT, const int NT, TMatrix& A, TMatrix& T )
{
	double ttime = omp_get_wtime();

	for (int tk=0; tk < min(MT,NT)-1; tk++ )
	{
		//////////////////////////////////////////////////////////////
		// Left
		//////////////////////////////////////////////////////////////
		{
			GEQRT( A(tk+1,tk), T(tk+1,tk) );

			#ifdef DEBUG
			cout << "GEQRT(" << tk+1 << "," << tk << "," << tk << ")\n";
			#endif
		}

		for (int tj=tk+1; tj < NT; tj++)
		{
			LARFB( PlasmaLeft, PlasmaTrans,   A(tk+1,tk), T(tk+1,tk), A(tk+1,tj) );

			#ifdef DEBUG
			#pragma omp critical
			cout << "LARFB_L(" << tk+1 << "," << tj << "," << tk << ")\n";
			#endif
		} // j-LOOP END

		for (int ti=tk+2; ti < MT; ti++)
		{
			{
				TSQRT( A(tk+1,tk), A(ti,tk), T(ti,tk) );

				#ifdef DEBUG
				cout << "TSQRT( A(" << tk+1 << "," << tk << "), A(" << ti << "," << tk << ") )\n";
				#endif
			}

			for (int tj=tk+1; tj < NT; tj++)
			{
				SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk+1,tj), A(ti,tj) );

				#ifdef DEBUG
				#pragma omp critical
				cout << "SSRFB_L( A(" << tk+1 << "," << tj << "), A(" << ti << "," << tj << ") )\n";
				#endif
			} // j-LOOP END
		} // i-LOOP END

		//////////////////////////////////////////////////////////////
		// Right
		//////////////////////////////////////////////////////////////
		{
			// Copy A(tk+1,tk)^T to A(tk,tk+1)
			for (int i=0; i<A(tk+1,tk)->m(); i++)
				for (int j=0; j<A(tk+1,tk)->n(); j++) {
					A(tk,tk+1)->Set_Val( j, i, A(tk+1,tk)->Get_Val(i,j) );
				}

			#ifdef DEBUG
			cout << "GEQRT(" << tk << "," << tk+1 << "," << tk << ")\n";
			#endif

		}

		for (int tj=tk+1; tj < NT; tj++)
		{
			LARFB( PlasmaRight, PlasmaNoTrans, A(tk+1,tk), T(tk+1,tk), A(tj,tk+1) );

			#ifdef DEBUG
			#pragma omp critical
			cout << "LARFB_R(" << tj << "," << tk+1 << "," << tk << ")\n";
			#endif
		} // j-LOOP END

		for (int ti=tk+2; ti < MT; ti++)
		{
			{
				// Copy A(tk+1,tk)^T to A(tk,tk+1)
				for (int i=0; i<A(tk+1,tk)->m(); i++)
					for (int j=0; j<A(tk+1,tk)->n(); j++) {
						A(tk,tk+1)->Set_Val( j, i, A(tk+1,tk)->Get_Val(i,j) );
					}

				// Copy A(ti,tk)^T to A(tk,ti)
				for (int i=0; i<A(ti,tk)->m(); i++)
					for (int j=0; j<A(ti,tk)->n(); j++) {
						A(tk,ti)->Set_Val( j, i, A(ti,tk)->Get_Val(i,j) );
					}

				#ifdef DEBUG
				cout << "TSQRT(" << tk << "," << ti << "," << tk << ")\n";
				#endif
			}

			for (int tj=tk+1; tj < NT; tj++)
			{
				SSRFB( PlasmaRight, PlasmaNoTrans, A(ti,tk), T(ti,tk), A(tj,tk+1), A(tj,ti) );

				#ifdef DEBUG
				#pragma omp critical
				cout << "SSRFB_R( A(" << tj << "," << tk+1 << "), A(" << tj << "," << ti << ") )\n";
				#endif
			} // j-LOOP END
		} // i-LOOP END

	} // k-LOOP END
}
