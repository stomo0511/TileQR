//
//  LeftLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <omp.h>
#include "CoreBlas.h"
#include "Progress.h"

#ifdef VTRACE
#include <vt_user.h>
#endif

#ifndef __Test__Min__
#define __Test__Min__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __Test__Min__

void tileQR( const int MT, const int NT, TMatrix< Tile<double> >& A, TMatrix< Tile<double> >& T )
{
  // Progress Table
  Progress_Table Pt( MT, NT, min(MT,NT) );
  Pt.Init();
  Pt.setIJK(0, 0, 0, NYET);
	
  ////////////////////////////////////////////////////////////////////////////
  // Left Looking tile QR
  #pragma omp parallel for schedule(static,1)
  for (int tk = 0; tk < NT; tk++) {
			
    for (int tl = 0; tl < min(MT,tk); tl++) {

      Pt.check_waitIJK(tl, tl, tl);	// Check for GEQRT_(tl,tl,tl)

      // LARFB
      {
        #ifdef VTRACE
	VT_TRACER("LARFB");
        #endif

	LARFB( PlasmaLeft, PlasmaTrans, 
	       A(tl,tl), T(tl,tl), A(tl,tk) );
      }

      for (int ti = tl+1; ti < MT; ti++) {

	Pt.check_waitIJK(ti, tl, tl);	// Check for TSQRT_(ti,tl,tl)

	// SSRFB
	{
          #ifdef VTRACE
	  VT_TRACER("SSRFB");
          #endif

	  SSRFB( PlasmaLeft, PlasmaTrans, 
		 A(ti,tl), T(ti,tl), A(tl,tk), A(ti,tk) );
	}
      } // End of i-loop
    } // End of l-loop
			
    if ( tk < MT ) {

      // GEQRT
      {				
        #ifdef VTRACE
	VT_TRACER("GEQRT");
        #endif

	GEQRT( A(tk,tk), T(tk,tk) );
      }			

      // Progress table update
      Pt.setIJK(tk, tk, tk, DONE);
				
      for (int ti = tk+1; ti < MT; ti++) {
	//
	// TSQRT
	//
	{					
          #ifdef VTRACE
	  VT_TRACER("TSQRT");
          #endif

	  TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );
	}

	// Progress table update
	Pt.setIJK(ti, tk, tk, DONE);
					
      } // End of i-loop
    } // End if ( k < MT )
  } // End of k-loop
  // Left Looking tile QR END
  ////////////////////////////////////////////////////////////////////////////
}
