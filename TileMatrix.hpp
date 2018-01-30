/*
 * TileMatrix.hpp
 *
 *  Created on: 2017/06/20
 *      Author: stomo
 */

#ifndef TILEMATRIX_HPP_
#define TILEMATRIX_HPP_

#include <cassert>

/*
 *  @class TileMatrix
*/
class TileMatrix
{
private:
	const int m_;		// number of lows of the matrix
	const int n_;		// number of columns of the matrix
	int mb_;			// number of lows of the tile
	int nb_;			// number of columns of the tile
	int mt_;			// number of low tiles
	int nt_;			// number of column tiles
	int ib_;			// inner block size

public:
	double *top_;			// pointer for top address

	/**
	 * Constructor
	 *
	 * @param m number of lows of the matrix
	 * @param n number of columns of the matrix
	 */
	TileMatrix(	const int m, const int n )
	: m_(m), n_(n), mb_(m), nb_(n), mt_(1), nt_(1), ib_(n)
	{
		assert( m > 0 && n > 0 );

		try
		{
			top_ = new double[ m * n ];
		}
		catch (char *eb)
		{
			std::cerr << "Can't allocate memory space for TileMatrix class: " << eb << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	/**
	 * Constructor
	 *
	 * @param m number of lows of the matrix
	 * @param n number of columns of the matrix
	 * @param mb number of lows of the tile
	 * @param nb number of columns of the tile
	 */
	TileMatrix(	const int m, const int n, const int mb, const int nb, const int ib )
	: m_(m), n_(n), mb_(m), nb_(n),
	  mt_(m_ % mb_ == 0 ? m_ / mb_ : m_ / mb_ + 1),
	  nt_(n_ % nb_ == 0 ? n_ / nb_ : n_ / nb_ + 1),
	  ib_(ib)
	{
		assert( m > 0 && n > 0 && mb > 0 && nb > 0 && ib > 0);
		assert( mb <= m && nb <= n );

		try
		{
			top_ = new double[ m * n ];
		}
		catch (char *eb)
		{
			std::cerr << "Can't allocate memory space for TileMatrix class: " << eb << std::endl;
			exit(EXIT_FAILURE);
		}
	}

	/**
	 * Destructor
	 */
	~TileMatrix()
	{
		delete [] top_;
	}

	/*
	 * Getters
	 */
	double* top() { return top_; }
	int m() const { return m_; }
	int n() const { return n_; }
	int mt() const { return mt_; }
	int nt() const { return nt_; }
	int ib() const { return ib_; }

	/*
	 * get mb of (ti,tj) tile
	 */
	int mb( const int ti, const int tj ) const
	{
		assert( ti < mt_ && tj < nt_ );

		return ( ti != mt_ -1 ) ? mb_ : (m_ % mb_);
	}

	/*
	 * get nb of (ti,tj) tile
	 */
	int nb( const int ti, const int tj ) const
	{
		assert( ti < mt_ && tj < nt_ );

		return ( tj != nt_ -1 ) ? nb_ : (n_ % nb_);
	}

	/*
	 *  Assign random numbers to the elements
	 *  @param seed seed of random number generator
	 */
	void Set_Rnd( const unsigned seed )
	{
		assert( seed >= 0 );

		srand(seed);
		for (int i = 0; i < m_ * n_; i++)
			top_[i] = (double)rand() / RAND_MAX;
	}

	// Set matrix to the identity matrix
/*	void Set_Iden()
	{
		for (int i=0; i<m_; i++)
			for (int j=0; j<n_; j++)
				top_[ i + j*m_ ] = ( i == j ) ? (double)(1) : (double)(0);
	}*/

	/*
	 * return pointer to the top address of tile (ti,tj)
	 *
	 * @param ti tile index
	 * @param tj tile index
	 *
	 * @return pointer to the top address of tile (ti,tj)
	 */
	double* ttop( const int ti, const int tj ) const
	{
		assert( ti >= 0 && ti < mt_ );
		assert( tj >= 0 && tj < nt_ );

		// column major x column major
		return top_ + ti* (mb_*nb_) + tj * (m_*nb_);
	}

/*	void makeI()
	{
		for(long int i=0; i<nol_; ++i){
			for(long int j=0; j<noc_; ++j){
				long int pos = 0;
				pos += (i/m_tsize_)*m_tsize_*noc_;    //タイルのi方向
				pos += (i%m_tsize_);              //タイル内のi方向
				if( i/m_tsize_ == mt_-1 ){
					pos += (nol_%m_tsize_ == 0 ? m_tsize_ : (nol_%m_tsize_))*n_tsize_*(j/n_tsize_);  //タイルのj方向
					pos += (nol_%m_tsize_ == 0 ? m_tsize_ : (nol_%m_tsize_))*(j%n_tsize_); //タイル内のj方向
				} else {
					pos += m_tsize_*n_tsize_*(j/n_tsize_);  //タイルのj方向
					pos += m_tsize_*(j%n_tsize_);
				}


				if( i==j ){
					top_[pos] = 1.0;
				}
				else
					top_[pos] = 0.0;
			}
		}
		return;
	}*/

};

#endif /* TILEMATRIX_HPP_ */
