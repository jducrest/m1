#include <math.h>

/** Multiply two matrices a(of size mxk) and b(of size kxn), and add the result to c(of size mxn)
 */
void multiplyMatrix(
    double *a,
    double *b,
    double *c,
    int m,
    int k,
    int n);

/** Create a matrix of size nrowsxncols.
 *  If init is set to "random", elements are initialized randomly. If it is "zero", elements are initialized to zero.
 */
void createMatrix(
    double **pmat,
    int nrows,
    int ncols,
    char *init);

/* Parallel SUMMA matrix-matrix multiplication of two matrices of size NxN
 */
void summa(int N)
{
    int rank, p;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &p);
    p = (int)sqrt(p);

    /* Create communicators */
    int row_color = rank / p;
    MPI_Comm row_comm;
    MPI_Comm_split(MPI_COMM_WORLD, row_color, rank, &row_comm);

    int col_color = rank % p + p;
    MPI_Comm col_comm;
    MPI_Comm_split(MPI_COMM_WORLD, col_color, rank, &col_comm);

    /* create matrices */
    double * Aloc, * Bloc, * Cloc;
    createMatrix(&Aloc, N / p, N / p, "random");
    createMatrix(&Bloc, N / p, N / p, "random");
    createMatrix(&Cloc, N / p, N / p, "zero");

    double * Atemp, * Btemp;
    createMatrix(&Atemp, N / p, N / p, "zero");
    createMatrix(&Btemp, N / p, N / p, "zero");

    int size = N * N / p / p;

    /* compute */
    int k;
    for (k = 0; k < p; ++k) {
        if (col_color == k + p)
            memcpy(Atemp, Aloc, size);

        if (row_color == k)
            memcpy(Btemp, Bloc, size);

        MPI_Bcast(Atemp, size, MPI_DOUBLE, k, row_comm);
        MPI_Bcast(Btemp, size, MPI_DOUBLE, k, col_comm);

        multiplyMatrix(Atemp, Btemp, Cloc, N / p, N / p, N / p);
    }

    /* cleanup */
    MPI_Comm_free(&row_comm);
    MPI_Comm_free(&col_comm);
    free(Aloc);
    free(Bloc);
    free(Cloc);
    free(Atemp);
    free(Btemp);
}

void shift(
    double *matrix,
    int shiftSize,
    char shiftDirection)
{
  int sendProcRank;
  int recvProcRank;
  MPI_Comm_rank(MPI_COM_WORLD,&rank);
  int i = rank/p;
  int j = rank%p;
  if (shiftDirection == 'l') { // Shift left; compute sendProcRank and recvProcRank accordingly
		sendProcRank = p*i + ((j-1)%p);
		recvProcRank = p*i + ((j+1)%p);
  } else if (shiftDirection == 'u') { // Shift up; compute sendProcRank and recvProcRank accordingly
		sendProcRank = p*((i-1)%p) + j;
		sendProcRank = p*((i+1)%p) + j;
  } else {
    printf("ERROR: Unknown shift direction %c.\n", shiftDirection);
  }
  // Now receive from recvProcRank while sending the matrix to sendProcRank
  MPI_ISend(matrix, size, MPI_DOUBLE, sendProcRank, 0, MPI_COMM_WORLD);
  MPI_Revc(matrix,size,MPI_DOUBLE,recvProcRank, 0, MPI_COMM_WORLD);
}

/* Parallel Cannon matrix-matrix multiplication of two matrices of size NxN
 */
void cannon(int N)
{
}
