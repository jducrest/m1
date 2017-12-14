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
	 qsort(E,M,sizeof(edge),sed_lex);
	n=0;k=0;
	while(n<N-1)
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
			n++; // we found a good edge... happy =)
			S[S[MAX(x,y)]]=MIN(x,y);
			printf("%d %d\n",E[k].i,E[k].j);
			//printf("%d %d %d\n",E[k].w,E[k].i,E[k].j);

		}
		k++;// we prepare to look into next edge
	}
