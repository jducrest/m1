/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algoName the name of the algorithm to be executed
 */

#include <limits.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
typedef struct 
{
int i;
int j;
int w;
} edge;

int weightcomp(const void* u,const void* v)
{
	if(((edge*)u)->w > ((edge*)v)->w)
		return 1;
	else
		return 0;
}

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
	 for(i=0;i<N;i++)
	 for(j=i+1;j<N;j++)
	 {
		if(adj[i*N+j] != 0 && adj[i*N+j]<D[j])
			D[j] = adj[i*N+j];
	 }

	 // main while
	 for(n=0;n<N-1;n++)
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
			if(T[i] && adj[i*N+imin]==D[imin])
				break;
		}
		// update T and D
		printf("%d %d\n",MIN(i,imin),MAX(i,imin));
		T[imin]=1;
		for(j=1;j<N;j++)
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
	 int i,j,k,n,min,imin;
	 // T is the tree, S is the representant set, E the set of edges
	 int* T=malloc(sizeof(int)*N);
	 int* S = malloc(sizeof(int)*N);
	 edge* E=malloc(sizeof(edge)*M);
	 // we initialize T and E
	 for(i=0;i<N;i++)
	 {
		T[i]=0;
		S[i]=i;
	 }
	 k = 0;
	 for(i=0;i<N;i++)
	 for(j=i;j<N;j++)
	 {
		if(adj[i*N+j]!=0)
		{
			E[k].i = i;
			E[k].j = j;
			E[k].w = adj[i*N+j];
			k++;
		}
	 }
	//printf("############## Edges before #################\n");
	 //for(j=0;j<M;j++)
		//printf("%d %d %d\n",E[j].i,E[j].j,E[j].w);
	 // we sort E
	 qsort(E,M,sizeof(edge),weightcomp);
	 //printf("################ Edges after sort #####################\n");
	 //for(j=0;j<M;j++)
	 //printf("%d %d %d\n",E[j].i,E[j].j,E[j].w);
 
	n=0;k=0;
	while(n<N-1 && k<M)
	{
		int x = E[k].i;
		int y = E[k].j;
		int w = E[k].w;
		k++;  // we prepare to look into next edge
		if(S[x]==S[y]) // loop or autoloop
			continue;
		n++; // we found a good edge... happy =)
		int m  = MIN(S[x],S[y]);
		for(i=0;i<N;i++) // we update the disjoint set structure "comme un sac"
			if(S[i]==S[x] || S[i]==S[y])
				S[i] = m;


		printf("%d %d\n",MIN(x,y),MAX(x,y));
	}

	 
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
