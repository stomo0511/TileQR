//
//  SuperM.cpp
//
//  Created by T. Suzuki on 2013/07/12.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <cassert>
#include <cstdlib>
#include "SuperM.h"

/**
 * Constructor
 *
 * @param m number of lows of the matrix
 * @param n number of columns of the matrix
 */
template <class _Tp>
SuperM<_Tp>::SuperM( const int m, const int n )
{
	assert( m > 0 );
	assert( n > 0 );
	
	m_ = m;
	n_ = n;
	top_ = NULL;
}

/**
 * Destructor
 */
template <class _Tp>
SuperM<_Tp>::~SuperM()
{
}

/**
 * Assign random numbers to the elements
 *
 * @param seed Seed of random number generator
 */
template <class _Tp>
void SuperM<_Tp>::Set_Rnd( const unsigned seed )
{
	assert( seed >= 0 );
	
	srand(seed);
	for (int i = 0; i < m_ * n_; i++)
		top_[i] = (_Tp)rand() / RAND_MAX;
}

/**
 * Assign random numbers to the elements ( for integer )
 *
 * @param seed Seed of random number generator
 */
template <>
void SuperM<int>::Set_Rnd( const unsigned seed )
{
	assert( seed >= 0 );
	
	srand(seed);
	for (int i = 0; i < m_ * n_; i++)
		top_[i] = rand();
}

/**
 * Set matrix to the identity matrix
 */
template <class _Tp>
void SuperM<_Tp>::Set_Iden()
{
	for (int i = 0; i < m_; i++)
		for (int j = 0; j < n_; j++)
			top_[ i + j * m_ ] = ( i == j ) ? (_Tp)(1) : (_Tp)(0);
}

/**
 * Set matrix to the identity matrix
 */
template <class _Tp>
void SuperM<_Tp>::Set_Zero()
{
	for (int i = 0; i < m_ * n_; i++)
		top_[ i ] = (_Tp)(0);
}

/**
 * Assign the value to (i,j) element
 *
 * @param i vertical index of the element
 * @param j horizontal index of the element
 * @param val element value
 */
template <class _Tp>
void SuperM<_Tp>::Set_Val( const int i, const int j, const _Tp val )
{
	assert( i >= 0 );	assert( i < m_ );
	assert( j >= 0 );	assert( j < n_ );
	
	top_[ i + j * m_ ] = val;
}

/**
 * Show elements to the standard output
 */
template <class _Tp>
void SuperM<_Tp>::Show_all() const
{
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n_; j++) {
			std::cout << top_[ i + j * m_ ] << " ";
		}
		std::cout << std::endl;
	}
	std::cout << std::endl;
}

/**
 * Operator overload =
 */
template <class _Tp>
SuperM<_Tp>& SuperM<_Tp>::operator=(const SuperM<_Tp>& T)
{
	assert( m_ == T.m_ );
	assert( n_ == T.n_ );
	
	for (int i = 0; i < m_ * n_; i++)
		top_[i] = T.top_[i];
	
	return *this;
}

/**
 * Operator overload []
 */
template <class _Tp>
_Tp& SuperM<_Tp>::operator[]( const int i ) const
{
	assert( i >= 0 );	assert( i < m_ * n_ );
	
	return top_[ i ];
}

template <class _Tp>
_Tp& SuperM<_Tp>::operator()( const int i, const int j ) const
{
	assert( i >= 0 );  assert( i < m_ );
	assert( j >= 0 );  assert( j < n_ );

	return top_[ i + j * m_ ];
}

/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 */
template <class _Tp>
void SuperM<_Tp>::File_Out( const char* fname )
{
	std::ofstream matf(fname);
	if (!matf) {
		std::cerr << "Unable to open " << fname << std::endl;
		exit(1);
	}
	
	matf << m_ << std::endl;
	matf << n_ << std::endl;
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n_; j++) {
			matf << top_[ i + j * m_ ] << " ";
		}
		matf << std::endl;
	}
	matf.close();
}



/**
 * Save matrix elements to the file
 *
 * @param fname data file name
 * @param dig number of output digit
 */
template <class _Tp>
void SuperM<_Tp>::File_Out( const char* fname, const unsigned dig )
{
	std::ofstream matf(fname);
	if (!matf) {
		std::cerr << "Unable to open " << fname << std::endl;
		exit(1);
	}
	
	matf << m_ << std::endl;
	matf << n_ << std::endl;
	matf.precision(dig);
	for (int i = 0; i < m_; i++) {
		for (int j = 0; j < n_; j++) {
			matf << top_[ i + j * m_ ] << " ";
		}
		matf << std::endl;
	}
	matf.close();
}

///**
// * Explicit instantiation of the template class
// */
template class SuperM<int>;
template class SuperM<float>;
template class SuperM<double>;
