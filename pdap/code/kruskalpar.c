/*********************************************/
	 /*                                           */
	 /*               KRUSKAL PAR                 */
    /*               goto kp                     */
	 /*********************************************/
	 
	 int Ns = procRank != numProcs - 1 ? ceil((float)N/numProcs) : N - ceil((float)N/numProcs)*(numProcs-1);
	 int offset = procRank * ceil((float)N / numProcs);

	 int i,j,k,n,min,imin;
	 // S is the representant set, E the set of edges
	 int*  S  = malloc(sizeof(int)*N);
	 edge* E  = (edge*) malloc(sizeof(edge)*MAX(M,N));
	 edge* E1 = (edge*) malloc(sizeof(edge)*N);
	 edge* E2 = (edge*) malloc(sizeof(edge)*N);
	 int first_time = 1;
	 int K;
	 n=0;
	 while(1)
	 {
	 if(first_time)
	 {
	 /* PHASE 1 : **********************************/
	 /* we initialize E using the adjacency matrix */
	 /**********************************************/
	 first_time = 0;
	 k = 0;
	 for(i=0;i<Ns;i++)
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
	 E[k].i = INT_MAX;
	 K = k;
	 }
	 else
	 {
	 /* PHASE 1 bis: ************************************************/
	 /* we initialize E using E1 and E2                             */
	 /* we do the return cards and see what happens method to get E */
	 /***************************************************************/
		i=0;j=0;k=0;
		while(1)
		{
			if(E1[i].i==INT_MAX)
			{
				while(E2[j].i!=INT_MAX)
				{
					E[k]=E2[j];j++;k++;
				}
				break;
			}
			if(E2[j].i==INT_MAX)
			{
				while(E1[i].i!=INT_MAX)
				{
					E[k] = E1[k];i++;k++;
				}
				break;
			}
			if(sed_lex(&E1[i],&E2[j]))
			{
				E[k]=E2[i];i++;
			}
			else
			{
				E[k]=E1[j];j++;
			}
			k++;
		}
		K=k;
	 }

	 /* PHASE 2 : **********************************/
	 /* we initialize S and go for it!             */
	 /**********************************************/

	 for(i=0;i<N;i++)
		S[i]=i;
	 qsort(E,K,sizeof(edge),sed_lex);
	i=0;k=0;
	while(k<K)
	{
		int x = E[k].i;
		int y = E[k].j;
		int w = E[k].w;
		while(S[x]!=S[S[x]])
			S[x] = S[S[x]];
		while(S[y]!=S[S[y]])
			S[y] = S[S[y]];
		if(S[x]!=S[y])
		{
			S[S[MAX(x,y)]]=MIN(x,y);
			E1[i]=E[k];
			i++; // we found a good edge... happy =)
			//printf("%d %d\n",MIN(x,y),MAX(x,y));
		}
		k++;// we prepare to look into next edge
	}
   	E1[i].i = INT_MAX; // this is a sign the table ended there
    if( ( 1 << n ) >= numProcs)
		break;
	if(procRank%(2<<n)==(1<<n)) // i'm the kind of guy who likes to send!
	{
		printf("I, %d try to send to %d!\n",procRank,procRank-(1<<n));
		MPI_Send(E1,N-1,edge_type,procRank-(1<<n),0,MPI_COMM_WORLD); // at most a tree so size O(N)
		printf("I, %d send my stuff!\n",procRank);
	}
	if(procRank%(2<<n)==0 && procRank+(1<<n) < numProcs) // i'm the kind of guy who likes to receive!
	{
		printf("I, %d try to recv from %d!\n",procRank,procRank+(1<<n));
		MPI_Recv(E2,N-1,edge_type,procRank+(1<<n),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
		printf("I, %d recv my stuff!\n",procRank);
	}
	MPI_Barrier(MPI_COMM_WORLD);
	n++;
	}

	if(procRank==0)
	{
		for(i=0;i<N-1;i++)
			printf("%d %d %d\n", E1[i].w, E1[i].i, E1[i].j);
	}
	
	free(S);
