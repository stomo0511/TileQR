//
//  LeftLooking
//
//  Created by T. Suzuki on 2022/12/23.
//

#include <cstdlib>
#include <omp.h>
#include "CoreBlas.h"

using namespace std;

#ifdef TRACE
extern void trace_cpu_start();
extern void trace_cpu_stop(const char *color);
extern void trace_label(const char *color, const char *label);
#endif

void tileQR( const int MT, const int NT, TMatrix< Tile<double> >& A, TMatrix< Tile<double> >& T )
{
      ////////////////////////////////////////////////////////////////////////////
      // Left Looking tile QR Task parallel version
      #pragma omp parallel
      {
            #pragma omp single
            {
                  for (int tk = 0; tk < NT; tk++)
                  {
                        for (int tl = 0; tl < min(MT,tk); tl++)
                        {
                              #pragma omp task depend(in:*A(tl,tl), *T(tl,tl)) depend(inout:*A(tl,tk))
                              {
                                    #ifdef TRACE
						trace_cpu_start();
						trace_label("Cyan", "LARFB");
						#endif

                                    LARFB( PlasmaLeft, PlasmaTrans, A(tl,tl), T(tl,tl), A(tl,tk) );

                                    #ifdef TRACE
						trace_cpu_stop("Cyan");
						#endif

                                    #ifdef DEBUG
                                    #pragma omp critical
                                    {
                                          cout << "LARFB(" << tl << "," << tk << ") : " << omp_get_thread_num() << "\n";
                                    }
                                    #endif
                              }

                              for (int ti = tl+1; ti < MT; ti++)
                              {
                                    #pragma omp task depend(in:*A(ti,tl), *T(ti,tl)) depend(inout:*A(tl,tk), *A(ti,tk))
                                    {
                                          #ifdef TRACE
							trace_cpu_start();
							trace_label("Blue", "SSRFB");
							#endif

                                          SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tl), T(ti,tl), A(tl,tk), A(ti,tk) );

                                          #ifdef TRACE
							trace_cpu_stop("Blue");
							#endif

                                          #ifdef DEBUG
                                          #pragma omp critical
                                          {
                                                cout << "SSRFB(" << ti << "," << tk << ") : " << omp_get_thread_num() << "\n";
                                          }
                                          #endif
                                    }
                              } // End of i-loop
                        } // End of l-loop
          
                        if ( tk < MT )
                        {
                              #pragma omp task depend(inout:*A(tk,tk)) depend(out:*T(tk,tk))
                              {
                                    #ifdef TRACE
					      trace_cpu_start();
					      trace_label("Red", "GEQRT");
					      #endif

                                    GEQRT( A(tk,tk), T(tk,tk) );

                                    #ifdef TRACE
					      trace_cpu_stop("Red");
					      #endif

                                    #ifdef DEBUG
                                    #pragma omp critical
                                    {
                                          cout << "GEQRT(" << tk << "," << tk << ") : " << omp_get_thread_num() << "\n";
                                    }
                                    #endif
                              }		

                              for (int ti = tk+1; ti < MT; ti++)
                              {
                                    #pragma omp task depend(inout:*A(tk,tk), *A(ti,tk)) depend(out:*T(ti,tk))
                                    {
                                          #ifdef TRACE
							trace_cpu_start();
							trace_label("Green", "TSQRT");
		 					#endif

                                          TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );

                                          #ifdef TRACE
							trace_cpu_stop("Green");
							#endif

                                          #ifdef DEBUG
                                          #pragma omp critical
                                          {
                                                cout << "TSQRT(" << ti << "," << tk << ") : " << omp_get_thread_num() << "\n";
                                          }
                                          #endif
                                    }
                              } // End of i-loop
                        } // End if ( k < MT )
                  } // End of k-loop
            } // End of single section
      } // End of parallel section
      // Left Looking tile QR END
      ////////////////////////////////////////////////////////////////////////////
}
