#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "mpi.h"

#include "cannon-solution.c"

void createMatrix(
    double **pmat,
    int nrows,
    int ncols,
    char *init)
{
  double *mat;
  int i;
  *pmat = mat = (double *)malloc(nrows * ncols * sizeof(mat[0]));
  if (strcmp(init, "random") == 0) { // Initialize the matrix elements randomly
    srand(time(0));
    for (i = 0; i < nrows * ncols; i++) { mat[i] = rand() / (double)RAND_MAX; }
  } else { // Set all matrix elements to zero
    memset(mat, 0, nrows * ncols * sizeof(mat[0]));
  }
}

/** Multiply two matrices a(of size mxk) and b(of size kxn), and add the result to c(of size mxn)
 */
void multiplyMatrix(
    double *a,
    double *b,
    double *c,
    int m,
    int k,
    int n)
{
  int im, ik, in;
  for (im = 0; im < m; im++) {
    for (in = 0; in < n; in++) {
      for (ik = 0; ik < k; ik++) {
        c[im + in * m] += a[im + ik * m] * b[ik + in * k];
      }
    }
  }
}

void printUsage()
{
  printf("Usage: mpirun -np [num-procs] ./cannon N [summa | cannon]\n");
}

int main(int argc, char **argv) {
  MPI_Init(&argc, &argv);

  int procRank;
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  if (argc < 3 && procRank == 0) {
    printUsage();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }
  MPI_Barrier(MPI_COMM_WORLD);

  int N = atoi(argv[1]);

  if (strcmp(argv[2], "cannon") == 0) { cannon(N); }
  else if (strcmp(argv[2], "summa") == 0) { summa(N); }
  else if (procRank == 0) {
    printf("Unknown algorithm name: %s.\n", argv[2]);
  }

  MPI_Finalize();
  return 0;
}
