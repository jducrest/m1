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
	 /*************************************/
	 /*                                   */
	 /*         PRIM'S ALGORITHM          */
	 /*                                   */
	 /*************************************/

	 int i,j,n;
	 int* T=malloc(sizeof(int)*N);
	 int* D=malloc(sizeof(int)*N);
	 int* Ne=malloc(sizeof(int)*N);
	 /*INIT**********************/
	 /* 0 is in the tree        */
	 /* everybody else is not   */
	 /* we update the distances */
	 /***************************/
	 T[0] = 1;
	 D[0] = INT_MAX;
	 Ne[0] = INT_MAX;
	 for(i=1;i<N;i++)
	 {
		T[i] = 0;
		if(adj[i+0*N] != 0)
		{
			D[i] = adj[i+0*N];
			Ne[i] = 0;
		}
		else
		{
			D[i] = INT_MAX;
			Ne[i] = INT_MAX;
		}
	 }
	 /* MAIN LOOP ******************************************/
	 /* find n-1 edges to build the tree                   */
	 /* we denote the edge by (u,v) where u is in the tree */
	 /* at the beginning, umin = 0, vmin = 0, dmin=INT_MAX */
	 /* this ensure that the first existing edge is picked */
	 /******************************************************/
	 for(n=0;n<N-1;n++)
	 {
		int umin = 0;  
		int vmin = 0;
		for(i=1;i<N;i++)
		{
			if( !T[i] && D[i]<D[vmin] ) // clear candidate
			{
				umin = Ne[i];
				vmin = i;
end_datatype		}
			if( !T[i] && D[i]==D[vmin] ) //potential candidate (need to check lex order)
			{
				int nei = Ne[i];
				int u = MIN(i,nei);
				int v = MAX(i,nei);
				if(u<MIN(umin,vmin) || (u==MIN(umin,vmin) && v<MAX(umin,vmin)))
				{
					umin = nei;
					vmin = i;
				}
			}
		}
		/* OUTPUT SOLUTION ****************/
		printf("%d %d\n",MIN(umin,vmin),MAX(umin,vmin));

		/* UPDATE *****************************************/
		/* we update T and D and V being careful with     */
		/* lex order                                      */
		/**************************************************/
		T[vmin]=1;
		for(j=1;j<N;j++)
		{
			if(adj[vmin*N+j]!=0 && adj[vmin*N+j]<D[j])
			{
				D[j] = adj[vmin*N+j];
				Ne[j] = vmin;
			}
			if(adj[vmin*N+j]!=0 && adj[vmin*N+j]==D[j])
			{
				Ne[j] = MIN(vmin,Ne[j]);
			}
		}
	 }

	 /* FREE **********************/
	 /* usual C boring stuff      */
	 /*****************************/
	 free(T);
	 free(D);
	 free(Ne);

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
	 qsort(E,M,sizeof(edge),weightcomp);
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

	 /*************************************/
	 /*                                   */
	 /*     PRIM'S PAR ALGORITHM          */
	 /*                                   */
	 /*************************************/
    if(! (N%p))
		printf("we are fucked !\n");
	 int sizeofdata = N*(N/p)
	 int Ns = N/p;
	 int* adjs = malloc(sizeof(int)*sizeofdata);
	 MPI_Scatter(adj,sizeofdata,MPI_INT,ajds,sizeofdata,MPI_INT,0,MPI_COMM_WORLD);

	 int i,j,n;
	 int* T=malloc(sizeof(int)*N);
	 int* D=malloc(sizeof(int)*Ns);
	 int* Ne=malloc(sizeof(int)*Ns);

// TO BE CONTINUED

	 /*INIT**********************/
	 /* 0 is in the tree        */
	 /* everybody else is not   */
	 /* we update the distances */
	 /***************************/
	 T[0] = 1;
	 D[0] = INT_MAX;
	 Ne[0] = INT_MAX;
	 for(i=1;i<N;i++)
	 {
		T[i] = 0;
		if(adj[i+0*N] != 0)
		{
			D[i] = adj[i+0*N];
			Ne[i] = 0;
		}
		else
		{
			D[i] = INT_MAX;
			Ne[i] = INT_MAX;
		}
	 }
	 /* MAIN LOOP ******************************************/
	 /* find n-1 edges to build the tree                   */
	 /* we denote the edge by (u,v) where u is in the tree */
	 /* at the beginning, umin = 0, vmin = 0, dmin=INT_MAX */
	 /* this ensure that the first existing edge is picked */
	 /******************************************************/
	 for(n=0;n<N-1;n++)
	 {
		int umin = 0;  
		int vmin = 0;
		for(i=1;i<N;i++)
		{
			if( !T[i] && D[i]<D[vmin] ) // clear candidate
			{
				umin = Ne[i];
				vmin = i;
			}
			if( !T[i] && D[i]==D[vmin] ) //potential candidate (need to check lex order)
			{
				int nei = Ne[i];
				int u = MIN(i,nei);
				int v = MAX(i,nei);
				if(u<MIN(umin,vmin) || (u==MIN(umin,vmin) && v<MAX(umin,vmin)))
				{
					umin = nei;
					vmin = i;
				}
			}
		}
		/* OUTPUT SOLUTION ****************/
		printf("%d %d\n",MIN(umin,vmin),MAX(umin,vmin));

		/* UPDATE *****************************************/
		/* we update T and D and V being careful with     */
		/* lex order                                      */
		/**************************************************/
		T[vmin]=1;
		for(j=1;j<N;j++)
		{
			if(adj[vmin*N+j]!=0 && adj[vmin*N+j]<D[j])
			{
				D[j] = adj[vmin*N+j];
				Ne[j] = vmin;
			}
			if(adj[vmin*N+j]!=0 && adj[vmin*N+j]==D[j])
			{
				Ne[j] = MIN(vmin,Ne[j]);
			}
		}
	 }

	 /* FREE **********************/
	 /* usual C boring stuff      */
	 /*****************************/
	 free(T);
	 free(D);
	 free(Ne);













  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
