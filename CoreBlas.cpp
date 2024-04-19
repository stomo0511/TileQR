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

/*
 * GEQRT conputes a QR factorization of a tile A: A = Q * R
 */
void GEQRT( const int mb, const int nb, const int ib, double *A, const int lda, double *T, const int ldt)
{
    assert( mb == nb );

    const int NB = max(lda,ldt);
	
    double* WORK = new double[ ib*NB ];
    double* TAU = new double[ NB ];
	
    plasma_core_dgeqrt( mb, mb, ib, A, lda, T, ldt, TAU, WORK );
	
    delete [] WORK;
    delete [] TAU;
}

/*
 * TSQRT conputes a QR factorization of a rectangular matrix formed by cupling (N x N) upper triangular tile A1 on top of (M x N) tile A2
 *
 */
void TSQRT( const int mb1, const int nb1, const int mb2, const int nb2, const int ib, double *A1, const int lda1, double *A2, const int lda2, double *T, const int ldt )
{
    assert( mb1 == nb1 );
    assert( nb1 == nb2 );
	
    const int NB = max(lda1,ldt);
	
    double* WORK = new double[ ib*NB ];
    double* TAU = new double[ NB ];
	
    plasma_core_dtsqrt( mb2, nb2, ib, A1, lda1, A2, lda2, T, ldt, TAU, WORK );
	
    delete [] WORK;
    delete [] TAU;
}

/*
 * LARFB updates (mb x nb) tile C with the transformation formed with A and T
 */
void LARFB( plasma_enum_t side, plasma_enum_t trans,
            const int mb, const int nb, const int kb, const int ib,
            double *A, const int lda, double *T, const int ldt, double *C, const int ldc )
{
    assert( (side==PlasmaLeft) || (side==PlasmaRight) );
    assert( (trans==PlasmaTrans) || (trans==PlasmaNoTrans) );
	
    // const int M = C->m();
    // const int N = C->n();
    // const int K = A->n();

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

/*
 * SSRFB updates (mb1 x nb1) tile C1 and (mb2 x nb2) tile C2 with the transformation formed with A and T
 */
void SSRFB( plasma_enum_t side, plasma_enum_t trans,
            const int mb1, const int nb1, const int mb2, const int nb2, const int kb, const int ib,
            double *A, const int lda, double *T, const int ldt, double *C1, const int ldc1, double *C2, const int ldc2 )
{
    assert( (side==PlasmaLeft) || (side==PlasmaRight) );
    assert( (trans==PlasmaTrans) || (trans==PlasmaNoTrans) );
	
    if (side == PlasmaRight)
        assert( mb2 == mb1);

    if (side == PlasmaLeft)
        assert( nb2 == nb1);

    // const int K = A->n();
	
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
                 WORK, LDWORK);
	
    delete [] WORK;
}

/*
 * dorgqr: genarates (M x N) orthogonal matrix Q: A = Q x R
 */
// void dorgqr( const TMatrix< Tile<double> > A,
//              const TMatrix< Tile<double> > T,
//              TMatrix< Tile<double> >& Q )
// {
//     assert( A.m() == Q.m() );
	
//     const int aMT = A.mt();
//     const int aNT = A.nt();
//     const int qMT = Q.mt();
//     const int qNT = Q.nt();

//     for (int tk = min(aMT, aNT)-1; tk+1 >= 1; tk--)
//     {
//         for (int ti = qMT - 1; ti > tk; ti--)
//         {
//             #pragma omp parallel for
//             for (int tj = tk; tj < qNT; tj++)
//                 SSRFB( PlasmaLeft, PlasmaNoTrans,
//                        A(ti,tk), T(ti,tk), Q(tk,tj), Q(ti,tj) );
//         }
//         #pragma omp parallel for
//         for (int tj = tk; tj < qNT; tj++)
//             LARFB( PlasmaLeft, PlasmaNoTrans,
//                    A(tk,tk), T(tk,tk), Q(tk,tj) );
//     }
// }

/*
 * TTQRT conputes a QR factorization of a rectangular matrix formed by cupling (N x N) upper triangular tile A1 on top of (M x N) upper trapezoidal tile A2
 */
// void TTQRT( Tile<double> *A1, Tile<double> *A2, Tile<double> *T )
// {
//     const int M = A2->m();
//     const int N = A1->n();
	
//     assert( N == A2->n() );
	
//     const int IB = A1->ib();
//     const int LDA1 = A1->m();
//     const int LDA2 = A2->m();
//     const int LDT = T->m();
	
//     const int NB = max(LDA1,LDT);
	
//     double* WORK = new double[ IB*NB ];
//     double* TAU = new double[ NB ];
	
//     plasma_core_dttqrt( M, N, IB,
//                  A1->top(), LDA1,
//                  A2->top(), LDA2,
//                  T->top(), LDT,
//                  TAU, WORK );
	
//     delete [] WORK;
//     delete [] TAU;
// }

