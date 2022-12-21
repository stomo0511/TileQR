#include <iostream>
#include <cstdlib>
#include <fstream>
#include <sys/time.h>
#include <lapacke.h>

using namespace std;

#define IDX2C(i,j,ld) (((j)*(ld))+(i))

// timer
double  my_clock()
{
  struct timeval tv;
  gettimeofday(&tv, NULL);
  return tv.tv_sec + (double)tv.tv_usec*1e-6;
}

// generate test matrix
void gen_rnd_mat( int m, int n, int lda, double* a )
{
  srand(20110609);
  for (int j=0; j<n; j++)
    for (int i=0; i<m; i++)
      a[IDX2C(i,j,lda)] = (double)(rand()) / RAND_MAX;
}

void file_out( int m, int n, int lda, double* a, char* fname )
{
  ofstream matf(fname);
  if (!matf) {
    cerr << "Unable to open " << fname << endl;
    exit(1);
  }

  matf.precision(15);
  for (int i=0; i<m; i++) {
    for (int j=0; j<n; j++)
      matf << a[IDX2C(i,j,m)] << " ";
    matf << endl;
  }
  matf.close();
}

// Test routine for house()
int main(int argc, char* argv[])
{
  if (argc < 4) {
    cerr << "Usage: LAPACK_QR [m] [n] [# of steps]\n";
    return EXIT_FAILURE;
  }

  int m = atoi(argv[1]);
  int n = atoi(argv[2]);
  int p = atoi(argv[3]);
  double time;

  if (m < n) {
    cerr << "Must be ( m >= n )\n";
    return EXIT_FAILURE;
  }

  double* a = new double[m*n];
  double* tau = new double[n];

  int info;

  // Generate Random Matrix
  gen_rnd_mat(m,n,m,a);
//   file_out(m,n,m,a,"a.dat");

  time = my_clock();
  for (int k=0; k<p; k++) {
    info = LAPACKE_dgeqrf( LAPACK_COL_MAJOR, m, n, a, m, tau );
  }
  time = my_clock() - time;
  cout << m << ", " << n << ", " << time/(double)(p) << endl;

//   file_out(m,n,m,a,"r.dat");

  delete[] tau;
  delete[] a;

  return 0;
}
