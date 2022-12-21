//
//  Tile.h
//
//  Created by T. Suzuki on 2013/07/08.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#ifndef __Tile__Tile__
#define __Tile__Tile__

#include "SuperM.h"

/**
 * Tile class
 */
template <class _Tp>
class Tile : public SuperM<_Tp> {

private:
	int ib_;   // inner blocking size
	
public:
	//  Constructor
	Tile( const int m, const int n, const int ib=1 );  // Allow Tile(m,n)
	Tile( const Tile& T );   // copy constructor

	//  Destructor
	~Tile();

	//  Show the inner blocking size
	int ib() const { return ib_; }
};

#endif /* defined(__Tile__Tile__) */
