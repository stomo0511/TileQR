//
//  Matrix.cpp
//
//  Created by T. Suzuki on 2013/07/10.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <cassert>
#include "Matrix.h"

/**
 * Constructor
 *
 * @param m number of lows of the matrix
 * @param n number of columns of the matrix
 */
template <class _Tp>
Matrix<_Tp>::Matrix( const int m, const int n ) : SuperM<_Tp>(m,n)
{
	this->top_ = new _Tp[ this->m_ * this->n_ ];

	for (int i = 0; i < this->m_ * this->n_; i++ )
		this->top_ [ i ] = (_Tp)(0);
}

/**
 * Copy Constructor
 *
 * @param T original object
 */
template <class _Tp>
Matrix<_Tp>::Matrix( const Matrix<_Tp>& T ) : SuperM<_Tp>(T)
{
	this->top_ = new _Tp[ this->m_ * this->n_ ];
	for (int i = 0; i < this->m_ * this->n_; i++)
		this->top_[i] = T.top_[i];
}

/**
 * Destructor
 */
template <class _Tp>
Matrix<_Tp>::~Matrix()
{
	delete[] this->top_;
}

/**
  * Explicit instantiation of the template class
  */
template class Matrix<int>;
template class Matrix<float>;
template class Matrix<double>;
