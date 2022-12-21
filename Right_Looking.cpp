//
//  RightLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include "CoreBlas.h"

#ifdef VTRACE
#include <vt_user.h>
#endif

#ifndef __Test__Min__
#define __Test__Min__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __Test__Min__

void tileQR( const int MT, const int NT, TMatrix< Tile<double> >& A, TMatrix< Tile<double> >& T )
{
  for (int tk=0; tk < min(MT,NT); tk++ ) {
    #pragma omp parallel
    {
      #pragma omp single
      {
        #ifdef VTRACE
	VT_TRACER("GEQRT");
        #endif

	GEQRT( A(tk,tk), T(tk,tk) );
      }
		
      #pragma omp for schedule(dynamic)
      for (int tj=tk+1; tj < NT; tj++) {
        #ifdef VTRACE
	VT_TRACER("LARFB");
        #endif

	LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );
      }

      for (int ti=tk+1; ti < MT; ti++) {
        #pragma omp single
	{
          #ifdef VTRACE
	  VT_TRACER("TSQRT");
          #endif

	  TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );
	}
			
        #pragma omp for schedule(dynamic)
	for (int tj=tk+1; tj < NT; tj++) {
          #ifdef VTRACE
	  VT_TRACER("SSRFB");
          #endif

	  SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );
	} // j-LOOP END
      } // i-LOOP END
    } // parallel section END
  } // k-LOOP END
}
