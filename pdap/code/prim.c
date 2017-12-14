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
			edge e1;
			edge e2;
			e1.i = Ne[i];
			e1.j = i;
			e1.w = D[i];
			e2.i=umin;
			e2.j=vmin;
			e2.w = D[vmin];
			if(!T[i] && sed_lex(&e2,&e1))
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
