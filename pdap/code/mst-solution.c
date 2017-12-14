/** Computing the Minimum Spanning Tree of a graph
 * @param N the number of vertices in the graph
 * @param M the number of edges in the graph
 * @param adj the adjacency matrix
 * @param algoName the name of the algorithm to be executed
 */

#include <limits.h>
#define MIN(a,b) (((a)<(b))?(a):(b))
#define MAX(a,b) (((a)>(b))?(a):(b))
#define W(u,v) adj[u*N+v]
typedef struct 
{
int i;
int j;
int w;
} edge;


int sed_lex(const void* a, const void* b)
{
	edge* a_ = (edge*)a;
	edge* b_ = (edge*)b;
	int u,u_,v,v_,w,w_;
	u  = MIN(a_->i,a_->j);
	v  = MAX(a_->i,a_->j);
	u_ = MIN(b_->i,b_->j);
	v_ = MAX(b_->i,b_->j);
	w  = a_->w;
	w_ = b_->w;

	if(w  >  w_)
		return 1;
	if(w == w_)
	{
		if(u > u_)
			return 1;
		if(u == u_ && v > v_)
			return 1;
	}
	return 0;
}
void edge_reduce(edge* in, edge* inout, int* len,MPI_Datatype* dptr)
{
	int i;
	int u,u_,v,v_,w,w_;
	u  = MIN(in[i].i,in[i].j);
	v  = MAX(in[i].i,in[i].j);
	u_ = MIN(inout[i].i,inout[i].j);
	v_ = MAX(inout[i].i,inout[i].j);
	w  = in[i].w;
	w_ = inout[i].w;

	for(i=0;i<*len;i++)
	{
	if(w  >  w_)
		continue;
	if(w == w_)
	{
		if(u > u_)
			continue;
		if(u == u_ && v > v_)
			continue;
	}
	inout[i] = in[i];
	}
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

  MPI_Datatype edge_type;
  MPI_Type_contiguous(3,MPI_INT,&edge_type);
  MPI_Type_commit(&edge_type);

  MPI_Op edge_red;
  MPI_Op_create( (MPI_User_function *) edge_reduce, 1, &edge_red );

  if (strcmp(algoName, "prim-seq") == 0) { // Sequential Prim's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
	  #include "prim.c"
  } else if (strcmp(algoName, "kruskal-seq") == 0) { // Sequential Kruskal's algorithm
    if (procRank == 0) {
      if (numProcs != 1) {
        printf("ERROR: Sequential algorithm is ran with %d MPI processes.\n", numProcs);
        MPI_Abort(MPI_COMM_WORLD, 1);
      }
    }
    #include "kruskal.c"

	 
  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

	 #include "primpar.c"

  } else if (strcmp(algoName, "kruskal-par") == 0) { // Parallel Kruskal's algorithm
    // BEGIN IMPLEMENTATION HERE
	 #include "kruskalpar.c"

  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
