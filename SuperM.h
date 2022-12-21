//
//  SuperM.h
//
//  Created by T. Suzuki on 2013/07/12.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#ifndef __Tile__SuperM__
#define __Tile__SuperM__

/**
 * Super Matrix class
 */
template <class _Tp>
class SuperM {

protected:
	_Tp* top_;	// pointer to matrix
	int m_;		// number of lows of the matrix (m) or lda
	int n_;		// number of colums of the matrix (n)

public:
	//  Constructor
	SuperM( const int m, const int n );

	//  Destructor
	virtual ~SuperM();

	//  Return the pointer
	_Tp* top() { return top_; }

	//  Show the matrix size
	int m() const { return m_; }
	int n() const { return n_; }

	// Assign random numbers to the elements
	void Set_Rnd( const unsigned seed );
	
	// Set matrix to the identity matrix
	void Set_Iden();
	
	// Set matrix to the zero matrix
	void Set_Zero();
	
	// Assign the value to (i,j) element
	void Set_Val( const int i, const int j, const _Tp val );
	
	// Show elements to the standard output
	void Show_all() const;

	// Operator overload
	SuperM& operator=( const SuperM& T );
	_Tp& operator[]( const int i ) const;
	_Tp& operator()( const int i, const int j ) const;
	
	// Save matrix elements to the file
	void File_Out( const char* fname );
	void File_Out( const char* fname, const unsigned dig );
};

#endif /* defined(__Tile__SuperM__) */
