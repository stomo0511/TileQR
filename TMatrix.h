//
//  TMatrix.h
//
//  Created by T. Suzuki on 2013/07/10.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#ifndef __Tile__TMatrix__
#define __Tile__TMatrix__

#include "Tile.h"

/**
 * TMatrix class
 */
template <class _Tp>
class TMatrix {
		
private:
	_Tp** top_;     // pointer to matrix pointer
	int m_;    // # lows of the matrix (m) or lda
	int n_;    // # colums of the matrix (n)
	int mb_;   // # lows of the tile
	int nb_;   // # columns of the tile
	int mt_;   // # low tiles
	int nt_;   // # column tiles
	
public:
	//  Constructor
	TMatrix( const int m, const int n, const int mb, const int nb, const int ib );
	TMatrix( const TMatrix& T );   // copy constructor
	
	//  Destructor
	~TMatrix();

	//  Show # of lows of the Matrix
	int m() const { return m_; }
	
	//  Show # of columns of the Matrix
	int n() const { return n_; }
	
	//  Show # of lows of the tile
	int mb() const { return mb_; }
	
	// Show # of columns of the tile
	int nb() const { return nb_; }

	//  Show # of tile lows
	int mt() const { return mt_; }

	//  Show # of tile columns
	int nt() const { return nt_; }
	
	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );

	// Set matrix to the identity matrix
	void Set_Iden();

	//  Operator overload
	_Tp* operator[]( const int i ) const;
	_Tp* operator()( const int i, const int j ) const;
	
	// Save matrix elements to the file
	void File_Out( const char* fname );
	void File_Out( const char* fname, const unsigned dig );

	// copy matrix elements to the array
	void Array_Copy( double *array );
};

#endif /* defined(__Tile__TMatrix__) */
