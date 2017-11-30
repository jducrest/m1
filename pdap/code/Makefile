all: mst

mst: mst-skeleton.c mst-solution.c
	mpicc -O3 mst-skeleton.c -lm -o mst

clean:
	rm mst

sp:
	mpirun -np 1 ./mst simple.in prim-seq

sk:
	mpirun -np 1 ./mst simple.in kruskal-seq
