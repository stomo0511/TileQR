//
//  DynamicScheduling
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//  Reviced in 2013/09/03/22:44
//

//#define COUT
#define TRAC
//#define ANIM

#include <iostream>
#include <queue>
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

	// Task Queue
	queue<Triplet> Qu;
	Triplet F(0,0,0);
	Qu.push(F);
	
	bool run_flag = true;
	// Definitions and Initialize　END
	//////////////////////////////////////////////////////////////////////
	
	double ttime = omp_get_wtime();
	//////////////////////////////////////////////////////////////////////
	// Dynamic Scheduling tile QR
	#pragma omp parallel private(F) firstprivate(ttime)
	{
		#ifdef TRAC
		double start_t;
		#endif

		bool my_flag = true;
		bool my_turn = false;
		bool input = false;
		unsigned char s;
		
		while (my_flag)
		{
			#pragma omp critical (Queue)
			{
				if (!Qu.empty())
				{
					F = Qu.front(); Qu.pop();  // Dequeue task;
					my_turn = true;
				}
				else
					my_turn = false;
			}

			if (my_turn)
			{
				int ti = F.getI();
				int tj = F.getJ();
				int tk = F.getK();
				
				if (tj == tk)
				{
					if (ti == tk)
					{
						// GEQRT
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

						// Enqueue TSQRT task
						if (ti+1 < MT)
						{
							input = false;

							#pragma omp critical (Progress)
							{
								s = Pt.getIJK( ti+1, tj, tk );
								s |= DIR_I;
								Pt.setIJK( ti+1, tj, tk, s );
		
								if (s == DONE)
									input = true;
							}

							if (input)
							{
								F.setIJK( ti+1, tj, tk );
		
								#pragma omp critical (Queue)
								Qu.push(F);
							}
						} // Enqueue TSQRT task End
						else {
							#pragma omp critical (End)
							run_flag = false;  // End flag <- true
						}
						
						// Enqueue LARFB task
						for (int j=tk+1; j<NT; j++)
						{
							input = false;

							#pragma omp critical (Progress)
							{
								s = Pt.getIJK( ti, j, tk );
								s |= DIR_J;
								Pt.setIJK( ti, j, tk, s );
		
								if (s == DONE)
									input = true;
							}

							if (input)
							{
								F.setIJK( ti, j, tk );
		
								#pragma omp critical (Queue)
								Qu.push(F);
							}
						} // Enqueue LARFB task End
					} // GEQRT END
					else
					{
						// TSQRT
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

						// Enqueue TSQRT task
						if (ti+1 < MT)
						{
							input = false;

							#pragma omp critical (Progress)
							{
								s = Pt.getIJK( ti+1, tj, tk );
								s |= DIR_I;
								Pt.setIJK( ti+1, tj, tk, s );
		
								if (s == DONE)
									input = true;
							}

							if (input)
							{
								F.setIJK( ti+1, tj, tk );
		
								#pragma omp critical (Queue)
								Qu.push(F);
							}
						} // Enqueue TSQRT task End
						else if (tj+1 >= NT)
						{
							#pragma omp critical (End)
							run_flag = false;  // End flag <- true
						}
	    
						// Enqueue SSRFB task
						for (int j=tk+1; j<NT; j++)
						{
							input = false;

							#pragma omp critical (Progress)
							{
								s = Pt.getIJK( ti, j, tk );
								s |= DIR_J;
								Pt.setIJK( ti, j, tk, s );
		
								if (s == DONE)
									input = true;
							}

							if (input)
							{
								F.setIJK( ti, j, tk );
		
								#pragma omp critical (Queue)
								Qu.push(F);
							}
						} // Enqueue SSRFB task End
					} // TSQRT END
				}
				else
				{
					if (ti == tk)
					{
						// LARFB
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

						// Enqueue SSRFB task
						if (ti+1 < MT)
						{
							input = false;

							#pragma omp critical (Progress)
							{
								s = Pt.getIJK( ti+1, tj, tk );
								s |= DIR_I;
								Pt.setIJK( ti+1, tj, tk, s );
		
								if (s == DONE)
									input = true;
							}

							if (input)
							{
								F.setIJK( ti+1, tj, tk );
		
								#pragma omp critical (Queue)
								Qu.push(F);
							}
						} // Enqueue SSRFB task End
					} // LARFB END
					else
					{
						// SSRFB
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

						// Enqueue GEQRT task
						input = false;

						#pragma omp critical (Progress)
						{
							s = Pt.getIJK( ti, tj, tk+1 );
							s |= DIR_K;
							Pt.setIJK( ti, tj, tk+1, s );
							
							if (s == DONE)
								input = true;
						}

						if (input)
						{
							F.setIJK( ti, tj, tk+1 );
	      
							#pragma omp critical (Queue)
							Qu.push(F);
						}
						// Enqueue GEQRT task End
						
						// Enqueue SSRFB task
						if (ti+1 < MT)
						{
							input = false;

							#pragma omp critical (Progress)
							{
								s = Pt.getIJK( ti+1, tj, tk );
								s |= DIR_I;
								Pt.setIJK( ti+1, tj, tk, s );
								
								if (s == DONE)
									input = true;
							}

							if (input)
							{
								F.setIJK( ti+1, tj, tk );
								
								#pragma omp critical (Queue)
								Qu.push(F);
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
