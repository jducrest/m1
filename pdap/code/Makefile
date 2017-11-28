all: mst

mst: mst-skeleton.c mst-solution.c
	mpicc -O3 mst-skeleton.c -lm -o mst

clean:
	rm mst

run:
	mpirun -np 1 ./mst graph.in prim-seq
