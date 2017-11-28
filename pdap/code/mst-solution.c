/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algoName the name of the algorithm to be executed
 */

#include <limits.h>

void computeMST(
    int N,
    int M,
    int *adj,
    char *algoName)
{
  int procRank, numProcs;
  MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  MPI_Comm_size(MPI_COMM_WORLD, &numProcs);

  if (strcmp(algoName, "prim-seq") == 0) { // Sequential Prim's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE
	 int i,j,n,min,imin;
	 int* T=malloc(sizeof(int)*N);
	 int* D=malloc(sizeof(int)*N);

	 //init
	 for(i=1;i<N;i++)
	 {
		T[i] = 0;
		D[i] = INT_MAX;
	 }

	 T[0] = 1;
	 D[0] = INT_MAX;
	 for(j=0;j<N;j++)
	 {
		if(adj[i*N+j]==1 && adj[i*N+j]<D[j])
			D[j] = adj[i*N+j];
	 }

	 // main while
	 for(n=0;n<N;n++)
	 {
		// find best vertex to add
		imin = 0;
		for(i=1;i<N;i++)
		{
			if( !T[i] && D[i]<D[imin] )
				imin = i;
		}
		// find relative edge
		for(i=0;i<N;i++)
		{
			if(adj[i*N+imin]==D[imin])
				break;
		}
		// update T and D
		printf("%d %d\n",i,imin);
		T[imin]=1;
		for(j=0;j<N;j++)
		{
			if(adj[imin*N+j]!=0 && adj[imin*N+j]<D[j])
				D[j] = adj[imin*N+j];
		}
	 }



	 free(T);
	 free(D);

  } else if (strcmp(algoName, "kruskal-seq") == 0) { // Sequential Kruskal's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
