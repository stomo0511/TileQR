/*
 * geqrf.cpp
 *
 *  Created on: 2016/07/21
 *      Author: stomo
 */


#include <iostream>
#include <omp.h>
#include <cassert>
#include <cstdlib>


#include "Check_Accuracy.hpp"

#include <lapacke.h>

using namespace std;

int main(int argc, const char * argv[])
{
	if (argc < 3)
	{
		cerr << "Usage: a.out [M] [N]\n";
		exit (1);
	}

	const int M =  atoi(argv[1]);  // n. of rows of the matrix
	const int N =  atoi(argv[2]);  // n. of columns of the matrix

	assert( M >= N );

	#ifdef DEBUG
	cout << "M = " << M << ", N = " << N << endl;
	cout << "clock precision = " << omp_get_wtick() << endl;
	#endif

	return EXIT_SUCCESS;
}
