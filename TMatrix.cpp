//
//  TMatrix.cpp
//
//  Created by T. Suzuki on 2013/07/10.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <cassert>
#include <iostream>
#include <fstream>
#include "Matrix.h"
#include "Tile.h"
#include "TMatrix.h"

/**
 * Constructor
 *
 * @param m number of low tiles of the TMtrix
 * @param n number of column tiles of the TMatrix
 */
template <> TMatrix< Tile<double> >::TMatrix( const int m, const int n, const int mb, const int nb, const int ib )
{
  assert( m > 0 );	assert( n > 0 );
  assert( mb > 0 );	assert( nb > 0 );
  assert( ib > 0 );
	
  m_ = m;
  n_ = n;
	
  mb_ = mb;
  nb_ = nb;

  mt_ = ( m_ % mb_ == 0 ) ? m_ / mb_ : m_ / mb_ + 1;
  nt_ = ( n_ % nb_ == 0 ) ? n_ / nb_ : n_ / nb_ + 1;
	
  top_ = new Tile<double>* [ mt_ * nt_ ];
	
  for ( int i = 0; i < mt_; i++ )
    for ( int j = 0; j < nt_; j++ ) {
      int tm = ( i != mt_ - 1 ) ? mb_ : m_ - i * mb_;
      int tn = ( j != nt_ - 1 ) ? nb_ : n_ - j * nb_;
      top_[ i + j * mt_ ] = new Tile<double> (tm,tn,ib);
    }
}

/**
 * Constructor
 *
 * @param m number of low tiles of the TMtrix
 * @param n number of column tiles of the TMatrix
 */
template <> TMatrix< Tile<double> >::TMatrix( const TMatrix< Tile<double> >& T )
{
  m_ = T.m_;
  n_ = T.n_;
  mb_ = T.mb_;
  nb_ = T.nb_;
  mt_ = T.mt_;
  nt_ = T.nt_;

  top_ = new Tile<double>* [ mt_ * nt_ ];

  for ( int i = 0; i < mt_; i++ )
    for ( int j = 0; j  < nt_; j++ )
      top_[ i + j * mt_ ] = T.top_[ i + j * mt_ ];
}

/**
 * Destructor
 */
template <> TMatrix< Tile<double> >::~TMatrix()
{
//	delete[] this->top_;
}

/**
 * Assign random numbers to the elements
 *
 * @param seed Seed of random number generator
 */
template <> void TMatrix< Tile<double> >::Set_Rnd( const unsigned seed )
{
  Matrix<double> Tmp(m_,n_);
	
  Tmp.Set_Rnd( seed );

  // (I,J) : Index of the elements of Matrix
  for (int I = 0; I < m_; I++) {
    for (int J = 0; J < n_; J++) {
      // (ti,tj) : Tile Index
      int ti = I / mb_;
      int tj = J / nb_;
      // (i,j) : Index of the elements of Tile
      int i = I % mb_;
      int j = J % nb_;
			
      top_[ ti + tj*mt_ ]->Set_Val(i, j, Tmp(I,J));
    }
  }
}

/**
 * Set matrix to the identity matrix
 */
template <> void TMatrix< Tile<double> >::Set_Iden()
{
  for (int i = 0; i < mt_; i++)
    for (int j = 0; j < nt_; j++)
      if (i == j)
	top_[ i + j * mt_ ]->Set_Iden();
      else
	top_[ i + j * mt_ ]->Set_Zero();
}

/**
 * Operator overload []
 */
template <> Tile<double>* TMatrix< Tile<double> >::operator[]( const int i ) const
{
  assert( i >= 0 );	assert( i < m_ * n_ );
	
  return top_[ i ];
}

/**
 * Operator overload ()
 */
template <> Tile<double>* TMatrix< Tile<double> >::operator()( const int i, const int j ) const
{
  assert( i >= 0 );  assert( i < m_ );
  assert( j >= 0 );  assert( j < n_ );
	
  return top_[ i + j * mt_ ];
}

/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 */
template <> void TMatrix< Tile<double> >::File_Out( const char* fname )
{
  Matrix<double> Tmp(m_,n_);
	
  // (I,J) : Index of the elements of Matrix
  for (int I = 0; I < m_; I++) {
    for (int J = 0; J < n_; J++) {
      // (ti,tj) : Tile Index
      int ti = I / mb_;
      int tj = J / nb_;
      // (i,j) : Index of the elements of Tile
      int i = I % mb_;
      int j = J % nb_;

      double val = top_[ ti + tj*mt_ ]->operator()(i, j);
      Tmp.Set_Val( I, J, val );
    }
  }
  Tmp.File_Out( fname );
}

/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 * @param dig number of output digit
 */
template <> void TMatrix< Tile<double> >::File_Out( const char* fname, const unsigned dig )
{
  Matrix<double> Tmp(m_,n_);
	
  // (I,J) : Index of the elements of Matrix
  for (int I = 0; I < m_; I++) {
    for (int J = 0; J < n_; J++) {
      // (ti,tj) : Tile Index
      int ti = I / mb_;
      int tj = J / nb_;
      // (i,j) : Index of the elements of Tile
      int i = I % mb_;
      int j = J % nb_;
			
      double val = top_[ ti + tj*mt_ ]->operator()(i, j);
      Tmp.Set_Val( I, J, val );
    }
  }
  Tmp.File_Out( fname, dig );
}

/**
 * copy matrix elements to the array
 *
 * @param array array
 */
template <> void TMatrix< Tile<double> >::Array_Copy( double *array )
{
  Matrix<double> Tmp(m_,n_);
	
  // (I,J) : Index of the elements of Matrix
  for (int I = 0; I < m_; I++) {
    for (int J = 0; J < n_; J++) {
      // (ti,tj) : Tile Index
      int ti = I / mb_;
      int tj = J / nb_;
      // (i,j) : Index of the elements of Tile
      int i = I % mb_;
      int j = J % nb_;

      array[ I + J*m_ ] = top_[ ti + tj*mt_ ]->operator()(i, j);
    }
  }
}
