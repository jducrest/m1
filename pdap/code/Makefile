all: mst

mst: mst-skeleton.c mst-solution.c
	mpicc -O3 mst-skeleton.c -lm -o mst

clean:
	rm mst

graphp:
	mpirun -np 1 ./mst graph.in prim-seq

graphk:
	mpirun -np 1 ./mst graph.in kruskal-seq

graphpp:
	mpirun -np 5 ./mst graph.in prim-par

graphkp:
	mpirun -np 5 ./mst graph.in kruskal-par

graphkp1:
	mpirun -np 1 ./mst graph.in kruskal-par

graphkp2:
	mpirun -np 2 ./mst graph.in kruskal-par
simplep:
	mpirun -np 1 ./mst simple.in prim-seq

simplek:
	mpirun -np 1 ./mst simple.in kruskal-seq

valoup:
	mpirun -np 1 ./mst valou.in prim-seq

valouk:
	mpirun -np 1 ./mst valou.in kruskal-seq

simplepp:
	mpirun -np 5 ./mst simple.in prim-par

valoupp:
	mpirun -np 10 ./mst valou.in prim-par
