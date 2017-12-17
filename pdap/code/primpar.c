	/*************************************/
	/*                                   */
	/*     PRIM'S PAR ALGORITHM          */
	/*     goto pp                       */
	/*************************************/
	int Ns = procRank != numProcs - 1 ? ceil((float)N/numProcs) : N - ceil((float)N/numProcs)*(numProcs-1);
	int Ns_gen = ceil((float)N/numProcs);
	int offset = procRank * ceil((float)N / numProcs);

	int i,j,n;
   edge e,emin;
	int* T=malloc(sizeof(int)*Ns);
	int* D=malloc(sizeof(int)*Ns);
	int* Ne=malloc(sizeof(int)*Ns);
		
	/*INIT**********************/
	/* 0 is in the tree        */
	/* everybody else is not   */
	/* we update the distances */
	/***************************/ 
	for(i=0;i<Ns;i++)
		for(j=0;j<N;j++)
			if(W(i,j)==0)
				W(i,j) = INT_MAX;
	
	W(0,offset) = INT_MAX; // this is used in case we can't fine an edge to add
	 
	// Now everybody initialize his small thingies !!! Ns is just a init value, index range from 0 to <Ns
	for(i=0;i<Ns;i++)
	{
		T[i] = 0;
		D[i] = W(i,0);
		Ne[i] = 0;
	}
	if(procRank==0)
		T[0] = 1;
	
	/* MAIN LOOP ******************************************/
	/* find n-1 edges to build the tree                   */
	/* we denote the edge by (u,v) where u is in the tree */
	/* at the beginning, umin = 0, vmin = 0, dmin=INT_MAX */
	/* this ensure that the first existing edge is picked */
	/******************************************************/
	 
	for(n=0;n<N-1;n++)
	{
		// we find a potential edge
		emin.i=INT_MAX;
		emin.j=INT_MAX;
		emin.w=INT_MAX;
		for(i=0;i<Ns;i++)
			if(!T[i])
			{	
				e.i = i+offset;	e.j = Ne[i];	e.w = W(i,Ne[i]);
			
				if(sed_lex(&emin,&e))
					emin=e;
			}
		
		MPI_Allreduce(&emin,&emin,1,edge_type,edge_red,MPI_COMM_WORLD);
		
		/* UPDATE *****************************************/
		/* we update T and D and V being careful with     */
		/* lex order                                      */
		/**************************************************/
		if(procRank==0)
			printf("%d %d\n",MIN(emin.i,emin.j),MAX(emin.i,emin.j));
		
		if(emin.i/Ns_gen==procRank) //it's somebody from this proc that was chosen
		   T[emin.i-offset]=1;

		for(j=0;j<Ns;j++)
		{
			if(W(j,emin.i) < D[j])
			{
				D[j] = W(j,emin.i);
				Ne[j] = emin.i;
			}
			if(W(j,emin.i)==D[j])
				Ne[j] = MIN(emin.i,Ne[j]);
		}
	 }

	 /* FREE **********************/
	 /* usual C boring stuff      */
	 /*****************************/
	 free(T);
	 free(D);
	 free(Ne);
