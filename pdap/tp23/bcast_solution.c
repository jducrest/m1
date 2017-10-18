// Already defined variables in the original source code. Do not uncomment any of those, just use them directly!
// Process rank 0 should be the source of the broadcast
//
// int num_procs;
// int rank;
// char *bcast_implementation_name:   the bcast implementation name (argument #1)
// int chunk_size:                    the chunk size (optional argument #2)
// int NUM_BYTES:                     the number of bytes to broadcast
// char *buffer:                      the buffer to broadcast
//
// The method names should be:
// default_bcast
// naive_bcast
// ring_bcast
// pipelined_ring_bcast
// asynchronous_pipelined_ring_bcast
// asynchronous_pipelined_bintree_bcast
//
// GOOD LUCK (gonna need it)!

if (strcmp(bcast_implementation_name, "default_bcast") == 0) 
{ 
	// Just calling the library routine.
	MPI_Bcast(buffer, NUM_BYTES, MPI_CHAR, 0, MPI_COMM_WORLD);
} 
else if (strcmp(bcast_implementation_name, "naive_bcast") == 0) 
{ 
	// Send to all processes one-by-one from the root.
	if(rank==0)
	{
		int i;
		for(i=1;i<num_procs;i++)
			MPI_Send(buffer,  NUM_BYTES, MPI_CHAR, i, 0, MPI_COMM_WORLD);
	} 
	else
	{
		MPI_Recv(buffer, NUM_BYTES, MPI_CHAR, 0, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
	}
}
else if (strcmp(bcast_implementation_name, "ring_bcast") == 0) {
	// P_i sends to P_i+1
	if(rank>0)
		MPI_Recv(buffer, NUM_BYTES, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

	if(rank<num_procs-1)
		MPI_Send(buffer,  NUM_BYTES, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD);
}
else if (strcmp(bcast_implementation_name, "pipelined_ring_bcast") == 0) {
	// P_i sends to P_i+1 but using pipelining
	int sent_bytes = 0;
	while(sent_bytes<NUM_BYTES)
	{
		if(sent_bytes+chunk_size>NUM_BYTES)
			chunk_size = NUM_BYTES-sent_bytes;

		if(rank>0)
			MPI_Recv(buffer+sent_bytes,chunk_size, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(rank<num_procs-1)
			MPI_Send(buffer+sent_bytes,chunk_size, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD);
		sent_bytes+=chunk_size;
	}
}
else if (strcmp(bcast_implementation_name, "asynchronous_pipelined_ring_bcast") == 0) {
	// P_i sends to P_i+1 but using pipelining and asynchronoucity
	int sent_bytes = 0;
	MPI_Request request[1000000];
	int i = 0;
	while(sent_bytes<NUM_BYTES)
	{
		if(sent_bytes+chunk_size>NUM_BYTES)
			chunk_size = NUM_BYTES-sent_bytes;

		if(rank>0)
			MPI_Recv(buffer+sent_bytes,chunk_size, MPI_CHAR, rank-1, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

		if(rank<num_procs-1)
		{
			MPI_Isend(buffer+sent_bytes,chunk_size, MPI_CHAR, rank+1, 0, MPI_COMM_WORLD,(request+i));
		}
		sent_bytes+=chunk_size;
		i++;
	}
	for(i=0;i<NUM_BYTES/chunk_size+1;i++)
		MPI_Wait((request+i),MPI_STATUS_IGNORE);
}
else if (strcmp(bcast_implementation_name, "asynchronous_pipelined_bintree_bcast") == 0) {
///////////////////////////////////////
//// NOT IMPLEMENTED //////////////////
// it should be theoretically better //
// ... but la flemme						 //
///////////////////////////////////////
}

