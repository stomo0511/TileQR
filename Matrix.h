//
//  Matrix.h
//
//  Created by T. Suzuki on 2013/07/08.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#ifndef __Tile__Matrix__
#define __Tile__Matrix__

#include "SuperM.h"

/**
 * Matrix class
 */
template <class _Tp>
class Matrix  : public SuperM<_Tp> {

public:
	//  Constructor
	Matrix( const int m, const int n );
	Matrix( const Matrix& T );   // copy constructor
	
	//  Destructor
	~Matrix();	
};

#endif /* defined(__Tile__Matrix__) */
