all: mst

mst: mst-skeleton.c mst-solution.c
	mpicc -O3 mst-skeleton.c -lm -o mst

clean:
	rm mst

simplep:
	mpirun -np 1 ./mst simple.in prim-seq

simplek:
	mpirun -np 1 ./mst simple.in kruskal-seq

valoup:
	mpirun -np 1 ./mst valou.in prim-seq

valouk:
	mpirun -np 1 ./mst valou.in kruskal-seq
