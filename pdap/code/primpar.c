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
		
		/* UPDATE *****************************************/
		/* we update T and D and V being careful with     */
		/* lex order                                      */
		/**************************************************/
		umin = e.i;
		vmin = e.j;
		if(procRank==0)
			printf("%d %d\n",MIN(umin,vmin),MAX(umin,vmin));
		
		if(vmin/Ns==procRank) //it's somebody from this proc that was chosen
	   {
		   T[vmin-offset]=1;
		}
		for(j=0;j<Ns;j++)
		{
			if(W(j,vmin)!=0 && W(j,vmin)<D[j])
			{
				D[j] = W(j,vmin);
				Ne[j] = vmin;
			}
			if(W(j,vmin)!=0 && W(j,vmin)==D[j])
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
