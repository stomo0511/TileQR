//
//  Progress.cpp
//  Tile
//
//  Created by T. Suzuki on 2013/08/18.
//  Copyright (c) 2013 T. Suzuki. All rights reserved.
//

#include <cassert>
#include "Progress.h"

#ifndef __Test__Min__
#define __Test__Min__

#define min(a,b) (((a)<(b)) ? (a) : (b))

#endif // __Test__Min__

using namespace std;

// Constructor
Progress_Table::Progress_Table( const int M, const int N, const int K )
{
  m_ = M;
  n_ = N;
  k_ = K;
  top_ = new unsigned char[ m_*n_*k_ ];
}


// Update the Table
void Progress_Table::setIJK( const int i, const int j, const int k, const unsigned char stat )
{
  assert( i>=0 );
  assert( j>=0 );
  assert( k>=0 );
	
  top_[ i + j*m_ + k*m_*n_ ] = stat;
}

// Get the Value of Table
unsigned char Progress_Table::getIJK( const int i, const int j, const int k )
{
  assert( i>=0 );
  assert( j>=0 );
  assert( k>=0 );
	
  return top_[ i + j*m_ + k*m_*n_ ];
}

// Check the Value of (i,j,k)
// Wait until it is set to DONE
void Progress_Table::check_waitIJK( const int i, const int j, const int k )
{
  while (1) {
    volatile bool gc;

    #pragma omp critical
    gc = top_[ i + j*m_ + k*m_*n_ ] == DONE ? true : false;

    if (gc)
      break;
  }
}

// Initialize the Table
void Progress_Table::Init()
{
  for (int k=0; k<min(m_,n_); k++)
    for (int i=k; i<m_; i++)
      for (int j=k; j<n_; j++) {
	unsigned char s;
	if (k==0)
	  s = DIR_K;
	else
	  s = NYET;
	if (i==k)
	  s |= DIR_I;
	if (j==k)
	  s |= DIR_J;

	setIJK( i, j, k, s );
      }
}

// Show the Table
void Progress_Table::Show()
{
  for (int k=0; k<min(m_,n_); k++) {
    cout << "k = " << k << endl;
    for (int i=0; i<m_; i++) {
      for (int j=0; j<n_; j++)
	cout << (int)(top_[ i + j*m_ + k*m_*n_ ]) << ", ";
      cout << endl;
    }
    cout << endl;
  }
}
