#include <plasma_core_blas.h>

#ifndef MKL
	#include <cblas.h>
	#include <lapacke.h>
#else
	#include <mkl.h>
#endif

void Check_Accuracy( const int M, const int N, double *mA, double *mQ, double *mR );
