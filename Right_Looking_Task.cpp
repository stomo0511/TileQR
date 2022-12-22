//
//  RightLooking Task Version
//
//  Created by T. Suzuki on 2022/12/21.
//

#include "CoreBlas.h"
#include <cstdlib>
#include <omp.h>

using namespace std;

#ifdef TRACE
extern void trace_cpu_start();
extern void trace_cpu_stop(const char *color);
extern void trace_label(const char *color, const char *label);
#endif

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
                    #pragma omp task depend(inout:Ap[tk][tk])
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
                }
            
                for (int tj=tk+1; tj < NT; tj++)
                {
                    #pragma omp task depend(in:Ap[tk][tk]) depend(inout:Ap[tk][tj])
                    {
                        #ifdef TRACE
						trace_cpu_start();
						trace_label("Cyan", "LARFB");
						#endif

                        LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );

                        #ifdef TRACE
						trace_cpu_stop("Cyan");
						#endif

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
                        #pragma omp task depend(inout:Ap[tk][tk]) depend(out:Ap[ti][tk])
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
                    }
                
                    for (int tj=tk+1; tj < NT; tj++)
                    {
                        #pragma omp task depend(in:Ap[ti][tk]) depend(inout:Ap[tk][tj],Ap[ti][tj])
                        {
                            #ifdef TRACE
							trace_cpu_start();
							trace_label("Blue", "SSRFB");
							#endif

                            SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );

                            #ifdef TRACE
							trace_cpu_stop("Blue");
							#endif

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
