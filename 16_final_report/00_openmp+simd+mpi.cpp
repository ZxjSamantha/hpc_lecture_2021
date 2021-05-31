#include "omp.h"
#include <mpi.h> // MPI (based on `module load intel-mpi`)
#include <cstdio>
#include <cstdlib> // drand48()
#include <cmath>
#include <vector>
#include <chrono> // time 
#include <immintrin.h>
using namespace std;
//bash 00_run.sh 
//1~ GFlops -> 20~ GFlops

int main(int argc, char** argv) {
  int size, rank;
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &size);
  MPI_Comm_rank(MPI_COMM_WORLD, &rank);

  //const int N = 256;
  //const int N = 512; 
  const int N = 1024;
  //const int N = 4096; 
  vector<float> A(N*N);
  vector<float> B(N*N);
  vector<float> C(N*N, 0);
  vector<float> subA(N*N/size); // size is the number of processes
  vector<float> subB(N*N/size);
  vector<float> subC(N*N/size, 0);
  
  //#pragma omp parallel for num_threads(4)
  #pragma omp parallel for collapse(2)
  for (int i=0; i<N; i++) {
    for (int j=0; j<N; j++) {
      A[N*i+j] = drand48();
      B[N*i+j] = drand48();
    }
  }


  int offset = N/size*rank; // begin = rank * (N/size) end = (rank + 1) * (N/size)
  
/* Parallel Section? */
#pragma omp parallel
{
  //int i,j; 
  #pragma omp sections 
  {
    for (int i=0; i<N/size; i++){
      for (int j=0; j<N; j++){
        subA[N*i+j] = A[N*(i+offset)+j];
      }
    }
  }
 
  #pragma omp sections
  {
    for (int i=0; i<N; i++){
      for (int j=0; j<N/size; j++){
        subB[N/size*i+j] = B[N*i+j+offset];
      }
    }
  }
}

  int recv_from = (rank + 1) % size; //receive 
  int send_to = (rank - 1 + size) % size; //send 

  double comp_time = 0, comm_time = 0;
  for(int irank=0; irank<size; irank++) {
    auto tic = chrono::steady_clock::now(); // record time 
    offset = N/size*((rank+irank) % size);
    #pragma omp parallel
    {
      int ic, jc, kc;
      #pragma omp for collapse(2)
      for (ic=0; ic<N/size; ic++){
        for (jc=0; jc<N/size; jc++){
          double temp = 0;
          /*reduction(+:temp) doesn't work well*/ 
          //#pragma omp parallel for reduction(+:temp)
          for (kc=0; kc<N; kc++){
            temp += subA[N*ic+kc] * subB[N/size*kc+jc];
            }  
          subC[N*ic+jc+offset] = temp; 
        }
      }
    }
    auto toc = chrono::steady_clock::now();
    comp_time += chrono::duration<double>(toc - tic).count();
    MPI_Send(&subB[0], N*N/size, MPI_FLOAT, send_to, 0, MPI_COMM_WORLD);
    MPI_Recv(&subB[0], N*N/size, MPI_FLOAT, recv_from, 0, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    tic = chrono::steady_clock::now();
    comm_time += chrono::duration<double>(tic - toc).count();
  }
  MPI_Allgather(&subC[0], N*N/size, MPI_FLOAT, &C[0], N*N/size, MPI_FLOAT, MPI_COMM_WORLD);
  
  omp_set_num_threads(4); // seems not matter? 
  
  #pragma omp parallel for //collapse(2)
  for (int i=0; i<N; i++){
    for (int j=0; j<N; j++){
      for (int k=0; k<N; k++){
        C[N*i+j] -= A[N*i+k] * B[N*k+j];
      }
    }
  }
  double err = 0;
  for (int i=0; i<N; i++)
    for (int j=0; j<N; j++)
      err += fabs(C[N*i+j]);
  if(rank==0) {
    double time = comp_time+comm_time;
    printf("N    : %d\n",N);
    printf("comp : %lf s\n", comp_time);
    printf("comm : %lf s\n", comm_time);
    printf("total: %lf s (%lf GFlops)\n",time,2.*N*N*N/time/1e9);
    printf("error: %lf\n",err/N/N);
  }
  MPI_Finalize();
}
