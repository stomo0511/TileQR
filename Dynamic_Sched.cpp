//
//  DynamicScheduling
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//  Reviced in 2013/09/03/22:44
//

#include <omp.h>
#include <queue>
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

  // Task Queue
  queue<Triplet> Qu;
  Triplet F(0,0,0);
  Qu.push(F);
	
  bool run_flag = true;
  // Definitions and Initialize　END
  //////////////////////////////////////////////////////////////////////
	
  //////////////////////////////////////////////////////////////////////
  // Dynamic Scheduling tile QR
  #pragma omp parallel private(F)
  {
    bool my_flag = true;
    bool my_turn = false;
    bool input = false;
    unsigned char s;
		
    while (my_flag) {

      {
	#ifdef VTRACE
	VT_TRACER("QUEUE");
	#endif

        #pragma omp critical (Queue)
	{
	  if (!Qu.empty()) {
	    F = Qu.front(); Qu.pop();  // Dequeue task;
	    my_turn = true;
	  }
	  else
	    my_turn = false;
	}
      }
			
      if (my_turn) {
	int ti = F.getI();
	int tj = F.getJ();
	int tk = F.getK();
				
	if (tj == tk) {
	  if (ti == tk) {
	    //
	    // GEQRT
	    //
	    {
              #ifdef VTRACE
	      VT_TRACER("GEQRT");
              #endif
	      GEQRT( A(tk,tk), T(tk,tk) );
	    }

	    // Enqueue TSQRT task
	    if (ti+1 < MT) {
	      input = false;

	      {
                #ifdef VTRACE
		VT_TRACER("PROGRESS");
                #endif

                #pragma omp critical (Progress)
		{
		  s = Pt.getIJK( ti+1, tj, tk );
		  s |= DIR_I;
		  Pt.setIJK( ti+1, tj, tk, s );
		
		  if (s == DONE)
		    input = true;
		}
	      }

	      {
                #ifdef VTRACE
		VT_TRACER("QUEUE");
                #endif

		if (input) {
		  F.setIJK( ti+1, tj, tk );
		
                  #pragma omp critical (Queue)
		  Qu.push(F);
		}
	      }
	    } // Enqueue TSQRT task End
	    else {
              #pragma omp critical (End)
	      run_flag = false;  // End flag <- true
	    }
						
	    // Enqueue LARFB task
	    for (int j=tk+1; j<NT; j++) {
	      input = false;

	      {
                #ifdef VTRACE
		VT_TRACER("PROGRESS");
                #endif

                #pragma omp critical (Progress)
		{
		  s = Pt.getIJK( ti, j, tk );
		  s |= DIR_J;
		  Pt.setIJK( ti, j, tk, s );
		
		  if (s == DONE)
		    input = true;
		}
	      }
	      {
                #ifdef VTRACE
		VT_TRACER("QUEUE");
                #endif

		{
                  #ifdef VTRACE
		  VT_TRACER("QUEUE");
                  #endif

		  if (input) {
		    F.setIJK( ti, j, tk );
		
                    #pragma omp critical (Queue)
		    Qu.push(F);
		  }
		}
	      }
	    } // Enqueue LARFB task End
	  } // GEQRT END
	  else {
	    //
	    // TSQRT
	    //
	    {
              #ifdef VTRACE
	      VT_TRACER("TSQRT");
              #endif
	      TSQRT( A(tk,tk), A(ti,tk), T(ti,tk) );
	    }

	    // Enqueue TSQRT task
	    if (ti+1 < MT) {
	      input = false;

	      {
                #ifdef VTRACE
		VT_TRACER("PROGRESS");
                #endif

                #pragma omp critical (Progress)
		{
		  s = Pt.getIJK( ti+1, tj, tk );
		  s |= DIR_I;
		  Pt.setIJK( ti+1, tj, tk, s );
		
		  if (s == DONE)
		    input = true;
		}
	      }
	      {
                #ifdef VTRACE
		VT_TRACER("QUEUE");
                #endif

		if (input) {
		  F.setIJK( ti+1, tj, tk );
		
                  #pragma omp critical (Queue)
		  Qu.push(F);
		}
	      }
	    } // Enqueue TSQRT task End
	    else if (tj+1 >= NT){

              #pragma omp critical (End)
	      run_flag = false;  // End flag <- true
	    }
	    
	    // Enqueue SSRFB task
	    for (int j=tk+1; j<NT; j++) {
	      input = false;

	      {
                #ifdef VTRACE
		VT_TRACER("PROGRESS");
                #endif

                #pragma omp critical (Progress)
		{
		  s = Pt.getIJK( ti, j, tk );
		  s |= DIR_J;
		  Pt.setIJK( ti, j, tk, s );
		
		  if (s == DONE)
		    input = true;
		}
	      }
	      {
                #ifdef VTRACE
		VT_TRACER("QUEUE");
                #endif

		if (input) {
		  F.setIJK( ti, j, tk );
		
                  #pragma omp critical (Queue)
		  Qu.push(F);
		}
	      } // Enqueue SSRFB task End
	    }
	  } // TSQRT END
	}
	else {
	  if (ti == tk) {
	    //
	    // LARFB
	    //
	    {
              #ifdef VTRACE
	      VT_TRACER("LARFB");
              #endif

	      LARFB( PlasmaLeft, PlasmaTrans, A(tk,tk), T(tk,tk), A(tk,tj) );
	    }
	    // Enqueue SSRFB task
	    if (ti+1 < MT) {
	      input = false;

	      {
                #ifdef VTRACE
		VT_TRACER("PROGRESS");
                #endif

                #pragma omp critical (Progress)
		{
		  s = Pt.getIJK( ti+1, tj, tk );
		  s |= DIR_I;
		  Pt.setIJK( ti+1, tj, tk, s );
		
		  if (s == DONE)
		    input = true;
		}
	      }

	      {
                #ifdef VTRACE
		VT_TRACER("QUEUE");
                #endif

		if (input) {
		  F.setIJK( ti+1, tj, tk );
		
                  #pragma omp critical (Queue)
		  Qu.push(F);
		}
	      }
	    } // Enqueue SSRFB task End
	  } // LARFB END
	  else {
	    //
	    // SSRFB
	    //
	    {
              #ifdef VTRACE
	      VT_TRACER("SSRFB");
              #endif

	      SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk,tj), A(ti,tj) );
	    }

	    // Enqueue GEQRT task
	    input = false;

	    {
              #ifdef VTRACE
	      VT_TRACER("PROGRESS");
              #endif

              #pragma omp critical (Progress)
	      {
		s = Pt.getIJK( ti, tj, tk+1 );
		s |= DIR_K;
		Pt.setIJK( ti, tj, tk+1, s );
							
		if (s == DONE)
		  input = true;
	      }
	    }
	    {
              #ifdef VTRACE
	      VT_TRACER("QUEUE");
              #endif

	      if (input) {
		F.setIJK( ti, tj, tk+1 );
	      
                #pragma omp critical (Queue)
		Qu.push(F);
	      }
	    }
	    // Enqueue GEQRT task End
						
	    // Enqueue SSRFB task
	    if (ti+1 < MT) {
	      input = false;

	      {
                #ifdef VTRACE
		VT_TRACER("PROGRESS");
                #endif

                #pragma omp critical (Progress)
		{
		  s = Pt.getIJK( ti+1, tj, tk );
		  s |= DIR_I;
		  Pt.setIJK( ti+1, tj, tk, s );
								
		  if (s == DONE)
		    input = true;
		}
	      }
	      {
                #ifdef VTRACE
		VT_TRACER("QUEUE");
                #endif

		if (input) {
		  F.setIJK( ti+1, tj, tk );
								
                  #pragma omp critical (Queue)
		  Qu.push(F);
		}
	      }
	    } // Enqueue SSRFB task End
	  } // SSRFB END
	}
      } // my_turn END
			
      #pragma omp critical (End)
      my_flag = run_flag;
    } // my_flag END
  } // End of outer most loop
  // Dynamic Scheduling QR END
  //////////////////////////////////////////////////////////////////////
}
