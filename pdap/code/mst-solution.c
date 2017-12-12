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

int sed_lex(int u, int v, int w, int u_, int v_, int w_) // answer first<second
{
	if(w<w_)
		return 1;
	if(w==w_ && MIN(u,v)<MIN(u_,v_))
		return 1;
	if(w==w_ && MIN(u,v)==MIN(u_,v_) && MAX(u,v)<MAX(u_,v_))
		return 1;

	return 0;
}


int weightcomp(const void* a, const void* b)
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
	 Ne[0] = INT_MAX;
	 D[0] = INT_MAX;
	 for(i=1;i<N;i++)
	 {
		T[i] = 0;
		if(W(0,i)!= 0)
		{
			D[i] = W(0,i);
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
		int u,v,w;
		for(i=1;i<N;i++)
		{
			if(!T[i] && sed_lex(Ne[i],i,D[i],umin,vmin,D[vmin]))
			{
				umin = Ne[i];
				vmin = i;
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
			if(W(vmin,j)!=0 && W(vmin,j)<D[j])
			{
				D[j] = W(vmin,j);
				Ne[j] = vmin;
			}
			if(W(vmin,j)!=0 && W(vmin,j)==D[j])
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
    /*****************************************/
	 /*                                       */
	 /*               KRUSKAL                 */
	 /*              goto k                   */
	 /*****************************************/
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
		if(W(i,j)!=0)
		{
			E[k].i = i;
			E[k].j = j;
			E[k].w = W(i,j);
			k++;
		}
	 }
	 qsort(E,M,sizeof(edge),weightcomp);
	 printf("QSORT:\n");
	 //for(i=0;i<M;i++)
	//	printf("%d %d %d\n", E[i].w, E[i].i, E[i].j);
	n=0;k=0;
	while(n<N-1)
	{
		int x = E[k].i;
		int y = E[k].j;
		int w = E[k].w;
		k++;// we prepare to look into next edge
		while(S[x]!=S[S[x]])
			S[x] = S[S[x]];
		while(S[y]!=S[S[y]])
			S[y] = S[S[y]];
		if(S[x]!=S[y])
		{
			n++; // we found a good edge... happy =)
			S[S[MAX(x,y)]]=MIN(x,y);
			printf("%d %d\n",MIN(x,y),MAX(x,y));

		}
	}

	 
  } else if (strcmp(algoName, "prim-par") == 0) { // Parallel Prim's algorithm
    // BEGIN IMPLEMENTATION HERE

	 /*************************************/
	 /*                                   */
	 /*     PRIM'S PAR ALGORITHM          */
	 /*     goto pp                       */
	 /*************************************/
	 int Ns = procRank != numProcs - 1 ? ceil((float)N/numProcs) : N - ceil((float)N/numProcs)*(numProcs-1);
	 int offset = procRank * ceil((float)N / numProcs);
	 // Now everybody got his adj small

	 int i,j,n;
    edge e;
	 int* T=malloc(sizeof(int)*Ns);
	 int* D=malloc(sizeof(int)*Ns);
	 int* Ne=malloc(sizeof(int)*Ns);
	 // thoses are the candidates that we send back to proc 0 for check
	 int* us = malloc(sizeof(int)*numProcs);
	 int* vs = malloc(sizeof(int)*numProcs);
	 int* ds = malloc(sizeof(int)*numProcs);
		
	 /*INIT**********************/
	 /* 0 is in the tree        */
	 /* everybody else is not   */
	 /* we update the distances */
	 /***************************/ 
	 W(0,offset) = INT_MAX; // this is used in case we can't fine an edge to add
	 
	 for(i=0;i<Ns;i++)
	 {
		T[i] = 0;
		D[i] = INT_MAX;
		Ne[i] = INT_MAX;
	 }
	 // Now everybody initialize his small thingies !!! Ns is just a init value, index range from 0 to <Ns
	 for(i=0;i<Ns;i++)
	 {
		if(W(i,0) != 0)
		{
			D[i] = W(i,0);
			Ne[i] = 0;
		}
	 }
	 if(procRank==0)
	 {
		T[0] = 1;
	 }
	 /* MAIN LOOP ******************************************/
	 /* find n-1 edges to build the tree                   */
	 /* we denote the edge by (u,v) where u is in the tree */
	 /* at the beginning, umin = 0, vmin = 0, dmin=INT_MAX */
	 /* this ensure that the first existing edge is picked */
	 /******************************************************/
	 
	 for(n=0;n<N-1;n++)
	 {
		int umin,vmin,dmin;
		// we find a potential edge
		for(i=0;i<Ns;i++)
		{
			if(!T[i])
			{
				umin = Ne[i];
				vmin = i+offset;
				break;
			}
		}
		if(i==Ns) // we couldn't find an edge give back error value
		{
			umin = INT_MAX;
			vmin = INT_MAX;
			dmin = INT_MAX;
			goto gatherpp;
		}
		// for the rest we compare to what we know
		for(i=i;i<Ns;i++)
		{
			if(!T[i] && D[i]<D[vmin-offset] ) // clear candidate
			{
				umin = Ne[i];
				vmin = i+offset;
			}
			if(!T[i] && D[i]==D[vmin-offset] ) //potential candidate (need to check lex order)
			{
				int nei = Ne[i];
				int u = MIN(i+offset,nei);
				int v = MAX(i+offset,nei);
				if(u<MIN(umin,vmin) || (u==MIN(umin,vmin) && v<MAX(umin,vmin)))
				{
					umin = nei;
					vmin = i+offset;
				}
			}
		}
		dmin = D[vmin-offset];
		gatherpp:
		
		e.i = umin;
		e.j = vmin;
		e.w = dmin;
		MPI_Allreduce(&e,&e,1,edge_type,edge_red,MPI_COMM_WORLD);
		/* noob method 
		// GATHER INFO 
		MPI_Gather(&umin,1,MPI_INT,us,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Gather(&vmin,1,MPI_INT,vs,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Gather(&dmin,1,MPI_INT,ds,1,MPI_INT,0,MPI_COMM_WORLD);
		int imin = 0;
		if(procRank==0)
		{
			for(i=1;i<numProcs;i++)
			{
				if(ds[i]<ds[imin])
					imin = i;
				if(ds[i] == ds[imin])
				if(MIN(us[i],vs[i]) < MIN(us[imin],vs[imin]) ||( MIN(us[i],vs[i]) == MIN(us[imin],vs[imin]) && MAX(us[i],vs[i])< MAX(us[imin],vs[imin])))
					imin = i;
			}
			// we found the best guy
			//  OUTPUT SOLUTION 
			
			if(procRank==0)
				printf("%d %d\n",MIN(us[imin],vs[imin]),MAX(us[imin],vs[imin]));
		}
		// we send it to everybody
		umin = us[imin];
		vmin = vs[imin];
		MPI_Bcast(&umin,1,MPI_INT,0,MPI_COMM_WORLD);
		MPI_Bcast(&vmin,1,MPI_INT,0,MPI_COMM_WORLD);
      //Now that everybody knows who it is we can print and update
      */
		/* UPDATE *****************************************/
		/* we update T and D and V being careful with     */
		/* lex order                                      */
		/**************************************************/
		umin = e.i;
		vmin = e.j;
		if(procRank==0)
			printf("%d %d\n",MIN(umin,vmin),MAX(umin,vmin));
			//printf("%d %d\n",MIN(us[imin],vs[imin]),MAX(us[imin],vs[imin]));
		
		if(vmin/Ns==procRank) //it's somebody from this proc that was chosen
	   {
		   T[vmin-offset]=1;
		}
		for(j=0;j<Ns;j++)
		{
			if(adj[vmin+j*N]!=0 && adj[vmin+j*N]<D[j])
			{
				D[j] = adj[vmin+j*N];
				Ne[j] = vmin;
			}
			if(adj[vmin+j*N]!=0 && adj[vmin+j*N]==D[j])
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
	 /*********************************************/
	 /*                                           */
	 /*               KRUSKAL PAR                 */
    /*               goto kp                     */
	 /*********************************************/

	 /* PHASE 1 : ******************************/
	 /* Do Kruskal on subset of edges          */
	 /******************************************/
	 
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
		if(W(i,j)!=0)
		{
			E[k].i = i;
			E[k].j = j;
			E[k].w = W(i,j);
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
		while(S[x]!=S[S[x]])
			S[x] = S[S[x]];
		while(S[y]!=S[S[y]])
			S[y] = S[S[y]];
		if(S[x]!=S[y])
		{
			n++; // we found a good edge... happy =)
			printf("%d %d\n",MIN(x,y),MAX(x,y));
		}
		int m = MIN(S[x],S[y]);
		for(i=0;i<N;i++)
			if(S[i]==S[x] || S[i] == S[y])
				S[i] = m;
	}


  } else { // Invalid algorithm name
    if (procRank == 0) {
      printf("ERROR: Invalid algorithm name: %s.\n", algoName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }
}
