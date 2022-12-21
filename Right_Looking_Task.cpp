//
//  RightLooking Task Version
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include "CoreBlas.h"
#include <cstdlib>
#include <omp.h>

#ifndef __Test__Min__
#define __Test__Min__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __Test__Min__

#ifdef VTRACE
#include <vt_user.h>
#endif

using namespace std;

void tileQR( const int MT, const int NT, TMatrix< Tile<double> >& A, TMatrix< Tile<double> >& T )
{
  //////////////////////////////////////////////////////////////////////
  // Progress Table
  int **Ap;

  Ap = (int **)malloc(sizeof(int *) * MT);
  for (int i=0; i<MT; i++)
      Ap[i] = (int *)malloc(sizeof(int) * NT);
  
  //////////////////////////////////////////////////////////////////////
  // Right Looking tile QR Task version
  #pragma omp parallel
  {
      #pragma omp single
      {
          for (int tk=0; tk < min(MT,NT); tk++ )
          {
              {
                  #ifdef VTRACE
                  VT_TRACER("GEQRT");
                  #endif

                  #pragma omp task depend(inout:Ap[tk][tk])
                  {
                      GEQRT( A(tk,tk), T(tk,tk) );

                      #ifdef DEBUG
                      #pragma omp critical
                      {
                          cout << "GEQRT(" << tk << "," << tk << ") : " << omp_get_thread_num() << "\n";
                      }
                      #endif
                  }
              }
		
              for (int tj=tk+1; tj < NT; tj++)
              {
                  #ifdef VTRACE
                  VT_TRACER("LARFB");
                  #endif

                  #pragma omp task depend(in:Ap[tk][tk]) depend(inout:Ap[tk][tj])
                  {
                      LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );

                      #ifdef DEBUG
                      #pragma omp critical
                      {
                          cout << "LARFB(" << tk << "," << tj << ") : " << omp_get_thread_num() << "\n";
                      }
                      #endif
                  }
              }

              for (int ti=tk+1; ti < MT; ti++)
              {
                  {
                      #ifdef VTRACE
                      VT_TRACER("TSQRT");
                      #endif

                      #pragma omp task depend(inout:Ap[tk][tk]) depend(out:Ap[ti][tk])
                      {
                          TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

                          #ifdef DEBUG
                          #pragma omp critical
                          {
                              cout << "TSQRT(" << ti << "," << tk << ") : " << omp_get_thread_num() << "\n";
                          }
                          #endif
                      }
                  }
			
                  for (int tj=tk+1; tj < NT; tj++)
                  {
                      #ifdef VTRACE
                      VT_TRACER("SSRFB");
                      #endif

                      #pragma omp task depend(in:Ap[ti][tk]) depend(inout:Ap[tk][tj],Ap[ti][tj])
                      {
                          SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );

                          #ifdef DEBUG
                          #pragma omp critical
                          {
                              cout << "SSRFB(" << ti << "," << tj << ") : " << omp_get_thread_num() << "\n";
                          }
                          #endif
                      }
                  } // j-LOOP END
              } // i-LOOP END
          } // k-LOOP END
      } // parallel section END
  }
  // Right Looking tile QR END
  //////////////////////////////////////////////////////////////////////
	
  for (int i=0; i<MT; i++)
      free(Ap[i]);
  free(Ap);
}
