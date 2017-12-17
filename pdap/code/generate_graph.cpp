#include <iostream>
#include <fstream>
#include <ctime>
#include <cmath>


int main(int argc, char **argv) {
    unsigned int N = std::atoi(argv[1]);
    int M = std::atoi(argv[2]);
    int W = std::atoi(argv[3]);
    std::ios_base::sync_with_stdio(false);
    unsigned int *edges = (unsigned int *)malloc(sizeof(unsigned int) * (N * (N-1) / 2));


    std::ofstream myfile;
  myfile.open (argv[4]);

    
    unsigned int index = 0;
    for (unsigned int i = 0; i < N; ++i) {
        for (unsigned int j = i+2; j < N; ++j) {
            edges[index++] = i*N+j;
        }
    }
    std::srand(std::time(nullptr)); 
    


    myfile<<N<<" "<<M<<"\n";
    for (int i = 0; i < N-1; ++i) {
        myfile<<i<<" "<<i+1<<" "<<std::rand() % W  + 1<<"\n";
        M--;
    }
    

    
    for (int i = 0; i < M; ++i) {
        unsigned int j = std::rand() % index;
        unsigned int a = edges[j] / N;
        unsigned int b = edges[j] % N;

        myfile<<a<<" "<<b<<" "<<std::rand() % W + 1<<"\n";
        edges[j] = edges[index];
        index--;
    }

    free(edges);
}
