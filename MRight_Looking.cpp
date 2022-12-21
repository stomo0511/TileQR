//
//  Modefied RightLooking
//
//  Created by T. Suzuki on 2015/07/24.
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

    for (int ti = 0; ti < MT; ti++)
    {
        #pragma omp parallel for schedule(dynamic)
        for (int tj = 0; tj < NT; tj++)
        {
            for (int tk = 0; tk <= min(ti,tj); tk++)
            {
                if (ti != tk)
                {
                    if ( tj == tk )
                    {
                        // TSQRT
                        {					
                            #ifdef VTRACE
                            VT_TRACER("TSQRT");
                            #endif

                            TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );
                        }
                        
                        // Progress table update
                        Pt.setIJK(ti, tk, tk, DONE);
                    }
                    else
                    {
                        // Check for TSQRT_(i,k,k)
                        Pt.check_waitIJK(ti, tk, tk);

                        // SSRFB
                        {
                            #ifdef VTRACE
                            VT_TRACER("SSRFB");
                            #endif

                            SSRFB( PlasmaLeft, PlasmaTrans, 
                                   A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );
                        }
                    }
                }
                else  // ti == tk
                {
                    if ( tj == tk )
                    {
                        // GEQRT
                        {				
                            #ifdef VTRACE
                            VT_TRACER("GEQRT");
                            #endif

                            GEQRT( A(tk,tk), T(tk,tk) );
                        }

                        // Progress table update
                        Pt.setIJK(tk, tk, tk, DONE);
                    }
                    else
                    {
                        // Check for TSQRT_(i,k,k)
                        Pt.check_waitIJK(tk, tk, tk);

                        // LARFB
                        {					
                            #ifdef VTRACE
                            VT_TRACER("LARFB");
                            #endif

                            LARFB( PlasmaLeft, PlasmaTrans, 
                                   A(tk,tk), T(tk,tk), A(tk,tj) );
                        }
                    }
                }
            }
        }
    }
    // Modefied Right Looking tile QR END
    ////////////////////////////////////////////////////////////////////////////
}
