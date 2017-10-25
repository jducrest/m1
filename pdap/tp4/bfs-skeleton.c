#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include <mpi.h>
#include <assert.h>

#include "bfs-solution.c"

void readGraph(
    char *fileName,
    int *pnVtx,
    int *pnEdge,
    int **padjBeg,
    int **padj)
{
  FILE *file = fopen(fileName, "r");
  if (file == NULL) {
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    if (rank == 0) {
      printf("ERROR: Unable to open the file %s.\n", fileName);
      MPI_Abort(MPI_COMM_WORLD, 1);
    }
  }

  int nVtx, nEdge;
  int *adjBeg, *adj;
  fscanf(file, " %d %d", pnVtx, pnEdge);
  nVtx = *pnVtx; nEdge = *pnEdge;
  int *edgeLeft, *edgeRight;
  edgeLeft = (int *) malloc(nEdge * sizeof(edgeLeft[0]));
  edgeRight = (int *) malloc(nEdge * sizeof(edgeRight[0]));
  *padjBeg = adjBeg = (int *) malloc((nVtx + 1) * sizeof(adjBeg[0]));
  memset(adjBeg, 0, (nVtx + 1) * sizeof(adjBeg[0]));
  *padj = adj = (int *) malloc(2 * nEdge * sizeof(adj[0]));
  
  int i;
  for (i = 0; i < nEdge; i++) {
    fscanf(file, " %d %d", edgeLeft + i, edgeRight + i);
    adjBeg[edgeLeft[i]]++;
    adjBeg[edgeRight[i]]++;
  }
  for (i = 1; i <= nVtx; i++) { adjBeg[i] += adjBeg[i - 1]; }
  for (i = 0; i < nEdge; i++) {
    adj[--adjBeg[edgeRight[i]]] = edgeLeft[i];
    adj[--adjBeg[edgeLeft[i]]] = edgeRight[i];
  }

  free(edgeLeft);
  free(edgeRight);
  fclose(file);
}

void printUsage() {
  int rank;
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  if (rank == 0) {
    printf("Usage: mpirun -np [num-procs] ./bfs [graph-file-name] [bfs-root-node]\n");
  }
}

void bfsSol(
    int N,
    int *adjBeg,
    int *adj,
    int *depth,
    int root)
{
  int visitedCount = 1;
  int prevVisitedCount = 0;
  while (prevVisitedCount < visitedCount) {
    prevVisitedCount = visitedCount;
    int curVtx;
    for (curVtx = 0; curVtx < N; curVtx++) {
      if (depth[curVtx] != -1) { continue; } // Vertex already visited
      int i;
      for (i = adjBeg[curVtx]; i < adjBeg[curVtx + 1]; i++) {
        if (depth[curVtx] == -1 || depth[adj[i]] + 1 < depth[curVtx]) { // Visit the neighbor, update the depth 
          depth[curVtx] = depth[adj[i]] + 1;
          visitedCount++;
        }
      }
    }
  }
}

int arrayEq(int N, int *arr1, int *arr2)
{
  int i;
  for (i = 0; i < N; i++) { if (arr1[i] != arr2[i]) { return 0; } }
  return 1;
}

int main(int argc, char **argv)
{
  MPI_Init(&argc, &argv);

  if (argc != 3) {
    printUsage();
    MPI_Abort(MPI_COMM_WORLD, 1);
  }

  int root = atoi(argv[2]);
  int procRank; MPI_Comm_rank(MPI_COMM_WORLD, &procRank);
  int nVtx, nEdge;
  int *adjBeg, *adj;
  int *depth, *depthSol;

  readGraph(argv[1], &nVtx, &nEdge, &adjBeg, &adj);

  depth = (int *) malloc(nVtx * sizeof(depth[0]));
  depthSol = (int *) malloc(nVtx * sizeof(depth[0]));
  memset(depth, -1, nVtx * sizeof(depth[0]));
  memset(depthSol, -1, nVtx * sizeof(depth[0]));
  depth[root] = depthSol[root] = 0;

  bfs(nVtx, adjBeg, adj, depth, root);
  bfsSol(nVtx, adjBeg, adj, depthSol, root);
  if (procRank == 0) {
    if (arrayEq(nVtx, depth, depthSol)) {
      printf("SUCCESS: BFS result is correct!\n");
    } else {
      printf("ERROR: BFS result is incorrect!\n");
    }
  }

  MPI_Finalize();
  return 0;
}
