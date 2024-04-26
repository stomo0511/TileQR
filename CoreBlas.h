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

void GEQRT( const int mb, const int nb, const int ib, double *A, const int lda, double *T, const int ldt);
void TSQRT( const int mb1, const int nb1, const int mb2, const int nb2, const int ib, double *A1, const int lda1, double *A2, const int lda2, double *T, const int ldt );
void LARFB( plasma_enum_t side, plasma_enum_t trans,
            const int mb, const int nb, const int kb, const int ib,
            double *A, const int lda, double *T, const int ldt, double *C, const int ldc );
void SSRFB( plasma_enum_t side, plasma_enum_t trans,
            const int mb1, const int nb1, const int mb2, const int nb2, const int kb, const int ib,
            double *A, const int lda, double *T, const int ldt, double *C1, const int ldc1, double *C2, const int ldc2 );

#endif /* defined(__Tile__CoreBlas__) */
