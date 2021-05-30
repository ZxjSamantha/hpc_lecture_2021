#include <cstdio>
#include <cstdlib>
#include <cmath>
#include "mpi.h"

//int main() {
int main(int argc, char** argv){
  const int N = 20;
  //const int np = 4;
  double x[N], y[N], m[N], fx[N], fy[N];
  MPI_Init(&argc, &argv); 
  int size, rank;
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);
  printf("rank: %d%d\n", rank, size);
  MPI_Finalize();
  int begin = rank * (N/size); 
  int end = (rank + 1) * (N/size);
  srand48(rank); 
  printf("%d %d %d\n", rank, begin, end); 

  for(int i=0; i<N; i++) {
    //if(i % (N / np) == 0) srand48(i / (N / np));
    x[i] = drand48();
    y[i] = drand48();
    m[i] = drand48();
    fx[i] = fy[i] = 0;
  }
  for(int i=0; i<N; i++) {
    for(int j=0; j<N; j++) {
      if(i != j) {
        double rx = x[i] - x[j];
        double ry = y[i] - y[j];
        double r = std::sqrt(rx * rx + ry * ry);
        fx[i] -= rx * m[j] / (r * r * r);
        fy[i] -= ry * m[j] / (r * r * r);
      }
    }
    printf("%d %g %g\n",i,fx[i],fy[i]);
  }
}
