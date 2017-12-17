all: mst smst

mst: mst-skeleton.c mst-solution.c
	mpicc -O3 mst-skeleton.c -lm -o mst

smst: mst-skeleton.c mst-solution.c
	smpicc -O3 mst-skeleton.c -lm -o smst

clean:
	rm mst;rm smst

densep:
	mpirun -np 1 ./mst dense.in prim-seq

densek:
	mpirun -np 1 ./mst dense.in kruskal-seq

densepp:
	mpirun -np 5 ./mst dense.in prim-par

densekp:
	mpirun -np 8 ./mst dense.in kruskal-par

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
	mpirun -np 32 ./mst valou.in prim-par


simplekp:
	mpirun -np 2 ./mst simple.in kruskal-par
valoukp:
	mpirun -np 10 ./mst valou.in kruskal-par
