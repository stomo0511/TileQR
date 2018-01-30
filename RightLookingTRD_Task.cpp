//
//  RightLooking
//
//  Created by T. Suzuki on 2014/01/05.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

//#define COUT
//#define DEBUG

#include <iostream>
#include <cstdlib>
#include <cassert>
#include <algorithm>

#include <CoreBlasTile.hpp>
#include <TMatrix.hpp>
#include <BMatrix.hpp>

using namespace std;

void tileCopy( BMatrix* A, BMatrix& B )
{
	assert( A->m() == B.m() );
	assert( A->m() == B.n() );
	assert( A->ib() == B.ib() );

	for (int i=0; i<A->m(); i++)
		for (int j=0; j<A->n(); j++)
			B.Set_Val(i,j,A->Get_Val(i,j));
}

void tileTransCopy( BMatrix* A, BMatrix& B )
{
	assert( A->m() == B.m() );
	assert( A->m() == B.n() );
	assert( A->ib() == B.ib() );

	for (int i=0; i<A->m(); i++)
		for (int j=0; j<A->n(); j++)
			B.Set_Val(j,i,A->Get_Val(i,j));
}

void tileTransCopy( BMatrix A, BMatrix* B )
{
	assert( A.m() == B->m() );
	assert( A.m() == B->n() );
	assert( A.ib() == B->ib() );

	for (int i=0; i<A.m(); i++)
		for (int j=0; j<A.n(); j++)
			B->Set_Val(j,i,A.Get_Val(i,j));
}

void tileTRD( const int MT, const int NT, TMatrix& A, TMatrix& T )
{
	const int NB = A(0,0)->m();
	const int IB = A(0,0)->ib();
	BMatrix tmp(NB,NB,IB), tmp1(NB,NB,IB);

	for (int tk=0; tk < min(MT,NT)-1; tk++ )
	{
		// (a)
		GEQRT( A(tk+1,tk), T(tk+1,tk) );
		#ifdef DEBUG
		cout << "GEQRT(" << tk+1 << "," << tk << "," << tk << ")\n";
		#endif

		// (b)
		LARFB( PlasmaLeft, PlasmaTrans,   A(tk+1,tk), T(tk+1,tk), A(tk+1,tk+1) );
		#ifdef DEBUG
		cout << "LARFB_L(" << tk+1 << "," << tk+1 << "," << tk << ")\n";
		#endif

		// (c)
		for (int ti=tk+1; ti < NT; ti++)
		{
			LARFB( PlasmaRight, PlasmaNoTrans, A(tk+1,tk), T(tk+1,tk), A(ti,tk+1) );
			#ifdef DEBUG
			cout << "LARFB_R(" << ti << "," << tk+1 << "," << tk << ")\n";
			#endif
		} // i-LOOP END

		//////////////////////////////////////////////////////////////
		for (int ti=tk+2; ti < MT; ti++)
		{
			// (d)
			TSQRT( A(tk+1,tk), A(ti,tk), T(ti,tk) );
			#ifdef DEBUG
			cout << "TSQRT( (" << tk+1 << "," << tk << "," << tk << "), (" << ti << "," << tk << "," << tk << ") )\n";
			#endif

			// (e) tj=tk+1
			tileTransCopy(A(ti,tk+1),tmp1);  // tmp1 <- A(ti,tk+1)^T

			SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), A(tk+1,tk+1), A(ti,tk+1) );
			#ifdef DEBUG
			cout << "SSRFB_L( (" << tk+1 << "," << tk+1 << "," << tk << "), (" << ti << "," << tk+1 << "," << tk << ") )\n";
			#endif

			// (f), (g)
			for (int tj=tk+2; tj<ti; tj++)
			{
				tileTransCopy(A(tj,tk+1),tmp);
				SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), &tmp, A(ti,tj) );
				tileTransCopy(tmp,A(tj,tk+1));
				#ifdef DEBUG
				cout << "SSRFB_L( (" << tj << "," << tk+1 << "," << tk << ")^T, (" << ti << "," << tj << "," << tk << ") )\n";
				#endif
			}

			// (h)
			SSRFB( PlasmaLeft, PlasmaTrans, A(ti,tk), T(ti,tk), &tmp1, A(ti,ti) );
			#ifdef DEBUG
			cout << "SSRFB_L( T(" << tk+1 << "," << ti << "," << tk << "), (" << ti << "," << ti << "," << tk << ") )\n";
			#endif

			// (i)
			SSRFB( PlasmaRight, PlasmaNoTrans, A(ti,tk), T(ti,tk), A(tk+1,tk+1), &tmp1 );
			#ifdef DEBUG
			cout << "SSRFB_R( (" << tk+1 << "," << tk+1 << "," << tk << "), T(" << tk+1 << "," << ti << "," << tk << ") )\n";
			#endif

			// (k)
			for (int tj=ti; tj < NT; tj++)
			{
				SSRFB( PlasmaRight, PlasmaNoTrans, A(ti,tk), T(ti,tk), A(tj,tk+1), A(tj,ti) );
				#ifdef DEBUG
				cout << "SSRFB_R( (" << tj << "," << tk+1 << "," << tk << "), (" << tj << "," << ti << "," << tk << ") )\n";
				#endif
			} // j-LOOP END
		} // i-LOOP END

		#ifdef DEBUG
		cout << endl;
		#endif
	} // k-LOOP END
}
