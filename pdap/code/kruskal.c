	/*****************************************/
	/*                                       */
	/*               KRUSKAL                 */
	/*              goto k                   */
	/*****************************************/
	int i,j,k,n;
	int* S = malloc(sizeof(int)*N);
	edge* E=malloc(sizeof(edge)*M);
	/*** INIT ********/
	/* init S and E  */
	/*****************/

	for(i=0;i<N;i++)
		S[i]=i;
	
	k = 0;
	for(i=0;i<N;i++)
	for(j=i+1;j<N;j++)
	{
		if(W(i,j)!=0)
		{
			E[k].i = i;
			E[k].j = j;
			E[k].w = W(i,j);
			k++;
		}
	}
	qsort(E,k,sizeof(edge),sed_lex);
	n=0;k=0;
	while(n<N-1)
	{
		int x = E[k].i;	int y = E[k].j;

		while(S[x]!=S[S[x]])
			S[x] = S[S[x]];
		while(S[y]!=S[S[y]])
			S[y] = S[S[y]];
		
		if(S[x]!=S[y])
		{
			S[S[MAX(x,y)]]=MIN(x,y);
			printf("%d %d \n",E[k].i,E[k].j);
			n++;
		}
		k++;
	}
	free(E);
	free(S);
