//
//  CoreBlas.cpp
//  Tile
//
//  Created by T. Suzuki on 2013/07/23.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <cassert>
#include <omp.h>
#include "CoreBlas.h"

using namespace std;

//
// GEQRT conputes a QR factorization of a tile A: A = Q * R
//
// @param[in] mb: row size of tile A
// @param[in] nb: col size of tile A
// @param[in] ib: inner block size
// @param[in,out] A: (mb x nb) tile matrix
// @param[in] lda: leading dimension of A
// @param[out] T: (ib x nb) upper triangular block reflector
// @param[in] ldt: leading dimension of T
//
void GEQRT( const int mb, const int nb, const int ib, 
    double *A, const int lda, double *T, const int ldt)
{
    assert( ib <= mb && ib <= nb  );

    const int NB = max(lda,ldt);
	
    double* WORK = new double[ ib*NB ];
    double* TAU = new double[ NB ];
	
    plasma_core_dgeqrt( mb, mb, ib, 
        A, lda, 
        T, ldt, 
        TAU, WORK );
	
    delete [] WORK;
    delete [] TAU;
}

//
// TSQRT conputes a QR factorization of a rectangular matrix formed by cupling (N x N) upper triangular tile A1 on top of (M x N) tile A2
// 
// @param[in] mb1: row size of tile A1
// @param[in] nb1: col size of tile A1
// @param[in] mb2: row size of tile A2
// @param[in] nb2: col size of tile A2
// @param[in] ib: inner block size
// @param[in,out] A1: (mb1 x nb1) tile matrix
// @param[in] lda1: leading dimension of A1
// @param[in,out] A2: (mb2 x nb2) tile matrix
// @param[in] lda2: leading dimension of A2
// @param[out] T: (ib x nb2) upper triangular block reflector
// @param[in] ldt: leading dimension of T
//
void TSQRT( const int mb1, const int nb1, const int mb2, const int nb2, const int ib, 
    double *A1, const int lda1, double *A2, const int lda2, double *T, const int ldt )
{
    assert( mb1 == nb1 );   // A1 should be square
    assert( nb1 == nb2 );   // A1 and A2 should have the same col size
	
    const int NB = max(lda1,ldt);
	
    double* WORK = new double[ ib*NB ];
    double* TAU = new double[ NB ];
	
    plasma_core_dtsqrt( mb2, nb2, ib, 
        A1, lda1, 
        A2, lda2, 
        T, ldt, TAU, WORK );
	
    delete [] WORK;
    delete [] TAU;
}

// 
// LARFB updates (mb x nb) tile C with the transformation formed with A and T
// 
// @param[in] side (PlasmaLeft or PlasmaRight)
// @param[in] trans (PlasmaTrans or PlasmaNoTrans)
// @param[in] mb: row size of tile C
// @param[in] nb: col size of tile C
// @param[in] kb: col size of tile A
// @param[in] ib: inner block size
// @param[in] A: (mb x kb) tile matrix
// @param[in] lda: leading dimension of A
// @param[in] T: (ib x kb) upper triangular block reflector
// @param[in] ldt: leading dimension of T
// @param[in,out] C: (mb x nb) tile matrix
// @param[in] ldc: leading dimension of C
//
void LARFB( plasma_enum_t side, plasma_enum_t trans,
            const int mb, const int nb, const int kb, const int ib,
            double *A, const int lda, double *T, const int ldt, double *C, const int ldc )
{
    assert( (side==PlasmaLeft) || (side==PlasmaRight) );
    assert( (trans==PlasmaTrans) || (trans==PlasmaNoTrans) );
	
    if (side == PlasmaLeft)
        assert( mb >= kb );
    else // (side == PlasmaRight)
        assert( nb >= kb );

    const int NB = max(lda,ldt);

    double* WORK = new double[ ib*NB ];

    plasma_core_dormqr( side, trans,
                 mb, nb, kb, ib,
                 A, lda,
                 T, ldt,
                 C, ldc,
                 WORK, nb );

    delete [] WORK;
}

//
// SSRFB updates (mb1 x nb1) tile C1 and (mb2 x nb2) tile C2 with the transformation formed with A and T
//
// @param[in] side (PlasmaLeft or PlasmaRight)
// @param[in] trans (PlasmaTrans or PlasmaNoTrans)
// @param[in] mb1: row size of tile C1
// @param[in] nb1: col size of tile C1
// @param[in] mb2: row size of tile C2
// @param[in] nb2: col size of tile C2
// @param[in] kb: col size of tile A
// @param[in] ib: inner block size
// @param[in] A: (mb1 x kb) tile matrix
// @param[in] lda: leading dimension of A
// @param[in] T: (ib x kb) upper triangular block reflector
// @param[in] ldt: leading dimension of T
// @param[in,out] C1: (mb1 x nb1) tile matrix
// @param[in] ldc1: leading dimension of C1
// @param[in,out] C2: (mb2 x nb2) tile matrix
// @param[in] ldc2: leading dimension of C2
//
void SSRFB( plasma_enum_t side, plasma_enum_t trans,
            const int mb1, const int nb1, const int mb2, const int nb2, const int kb, const int ib,
            double *A, const int lda, double *T, const int ldt, double *C1, const int ldc1, double *C2, const int ldc2 )
{
    assert( (side==PlasmaLeft) || (side==PlasmaRight) );
    assert( (trans==PlasmaTrans) || (trans==PlasmaNoTrans) );
	
    if (side == PlasmaRight)
        assert( mb2 == mb1);
    else // side == PlasmaLeft
        assert( nb2 == nb1);

    int LDWORK;
    if (side == PlasmaLeft)
        LDWORK = ib;
    else // side == PlasmaRight
        LDWORK = mb1;
	
    int WSIZE;
    if (side == PlasmaLeft)
        WSIZE = nb1;
    else // side == PlasmaRight
        WSIZE = ib;
	
    double* WORK = new double[ LDWORK * WSIZE ];

    plasma_core_dtsmqr( side, trans,
                 mb1, nb1, mb2, nb2, kb, ib,
                 C1, ldc1,
                 C2, ldc2,
                 A, lda,
                 T, ldt,
                 WORK, LDWORK );
	
    delete [] WORK;
}
