//
//  CoreBlas.h
//  Tile
//
//  Created by T. Suzuki on 2013/07/23.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#ifndef __Tile__CoreBlas__
#define __Tile__CoreBlas__

#include <iostream>
#include <plasma.h>
#include <plasma_core_blas.h>

#include "TMatrix.h"

void GEQRT( Tile<double> *A, Tile<double> *T );
void TSQRT( Tile<double> *A1, Tile<double> *A2, Tile<double> *T );
void TTQRT( Tile<double> *A1, Tile<double> *A2, Tile<double> *T );
void LARFB( plasma_enum_t side, plasma_enum_t trans,
			Tile<double> *A, Tile<double> *T, Tile<double> *C );
void SSRFB( plasma_enum_t side, plasma_enum_t trans,
		   Tile<double> *A, Tile<double> *T, Tile<double> *C1, Tile<double> *C2 );
void STRFB( plasma_enum_t side, plasma_enum_t trans,
		   Tile<double> *A, Tile<double> *T, Tile<double> *C1, Tile<double> *C2 );

void dorgqr( const TMatrix< Tile<double> > A, const TMatrix< Tile<double> > T, TMatrix< Tile<double> >& Q );

#endif /* defined(__Tile__CoreBlas__) */
