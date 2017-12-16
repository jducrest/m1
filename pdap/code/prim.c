	/*************************************/
	/*                                   */
	/*         PRIM'S ALGORITHM          */
	/*                                   */
	/*************************************/

	int i,j,n,umin,vmin;
	edge e,emin;
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
	/**************************************************/
	/* way easier to see it as a graph with           */
	/* infinite distance between non-connected edges. */
	/**************************************************/
	for(i=0;i<N;i++)
		for(j=0;j<N;j++)
			if(W(i,j)==0)
				W(i,j)=INT_MAX;

	for(i=1;i<N;i++)
	{
		T[i] = 0;
		D[i] = W(0,i);
		Ne[i] = 0;
	}
	/* MAIN LOOP ******************************************/
	/* find n-1 edges to build the tree                   */
	/* we denote the edge by (u,v) where u is in the tree */
	/* at the beginning, umin = 0, vmin = 0, dmin=INT_MAX */
	/* this ensure that the first existing edge is picked */
	/******************************************************/
	for(n=0;n<N-1;n++)
	{
		emin.w=INT_MAX;
		for(i=1;i<N;i++)
		{
			e.i = Ne[i];
			e.j = i;
			e.w = D[i];
			if(!T[i] && sed_lex(&emin,&e))
				emin = e;
		}
		/* OUTPUT SOLUTION ****************/
		printf("%d %d\n",MIN(emin.i,emin.j),MAX(emin.i,emin.j));

		/* UPDATE *****************************************/
		/* we update T and D and V being careful with     */
		/* lex order                                      */
		/**************************************************/
		i = emin.j;
		T[i]=1;
		for(j=1;j<N;j++)
		{
			if(W(i,j)<D[j])
			{
				D[j] = W(i,j);
				Ne[j] = i;
			}
			if(W(i,j)==D[j])
			{
				Ne[j] = MIN(i,Ne[j]);
			}
		}
	 }

	 /* FREE **********************/
	 /* usual C boring stuff      */
	 /*****************************/
	 free(T);
	 free(D);
	 free(Ne);
