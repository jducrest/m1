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
	edge* E  = (edge*) malloc(sizeof(edge)*MAX(M,N+1));
	edge* E1 = (edge*) malloc(sizeof(edge)*(N+1));
	edge* E2 = (edge*) malloc(sizeof(edge)*(N+1));
	int first_time = 1;
	int K;
	n=0;
	E1[0].i = INT_MAX;
	E2[0].i = INT_MAX;
	/* PHASE 1 : **********************************/
	/* we initialize E using the adjacency matrix */
	/**********************************************/
	first_time = 0;
	k = 0;
	for(i=0;i<Ns;i++)
	for(j=offset+i;j<N;j++)
	{
		if(W(i,j)!=0)
		{
			E[k].i = MIN(i+offset,j);
			E[k].j = MAX(i+offset,j);
			E[k].w = W(i,j);
			k++;
		}
	}
	E[k].i = INT_MAX;
	K = k;
	goto skipphaseone;
	while(1)
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
					E[k] = E1[i];i++;k++;
				}
				break;
			}
			if(sed_lex(&E1[i],&E2[j]) == sed_lex(&E2[j],&E1[i]))
			{
				E[k] = E2[j];i++;j++;	
			}
			else if(sed_lex(&E1[i],&E2[j])==1)
			{
				E[k]=E2[j];j++;
			}
			else
			{
				E[k]=E1[i];i++;
			}
			k++;
		}
		K=k;

		skipphaseone:
		/* PHASE 2 : **********************************/
		/* we initialize S and go for it!             */
		/**********************************************/
		if(procRank%(1<<n)==0)
		{
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
			}
			k++;// we prepare to look into next edge
		}
   	E1[i].i = INT_MAX; // this is a sign the table ended there
   }
	if( ( 1 << n ) >= numProcs)
		break;
	if(procRank%(2<<n)==(1<<n)) // i'm the kind of guy who likes to send!
	{
		MPI_Send(E1,N,edge_type,procRank-(1<<n),0,MPI_COMM_WORLD); // at most a tree so size O(N)
		break;
	}
	if(procRank%(2<<n)==0 && procRank+(1<<n) < numProcs) // i'm the kind of guy who likes to receive!
	{
		MPI_Recv(E2,N,edge_type,procRank+(1<<n),0,MPI_COMM_WORLD,MPI_STATUS_IGNORE);
	}
	n++;
	}
	if(procRank==0)
	{
		for(i=0;i<N-1;i++)
			printf("%d %d\n",E1[i].i,E1[i].j);
	}
	
	free(S);
