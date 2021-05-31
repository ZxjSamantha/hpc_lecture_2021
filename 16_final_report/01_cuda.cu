#include <mpi.h>
#include <cmath>
#include <cstdlib>
#include <cstdio>
#include <chrono>
using namespace std;
//bash 01_run_cuda.sh
//N = 256, 80~ GFlops
//N = 512, 210~ GFlops
//N = 1024, 350~ GFlops

#define M 256 
//number of rows in the matrix 
// N >= M

__global__ void matmul(float *A, float *B, float *C, int N){
    int i = blockIdx.y ; 
    int j = threadIdx.x + blockDim.x * blockIdx.x;
    float sum = 0;
    for (int k=0; k<N; k++){
        sum += A[N*i+k] * B[N*k+j]; 
    }
    C[N*i+j] = sum; 
}

int main(int argc, char **argv){
    int mpisize, mpirank;
    MPI_Init(&argc, &argv);
    MPI_Comm_size(MPI_COMM_WORLD, &mpisize);
    MPI_Comm_rank(MPI_COMM_WORLD, &mpirank);
    
    //const int N = 512;
    const int N = 1024; 
    int matsize = N * N * sizeof(float);
    float *A, *B, *C;   
    cudaMallocManaged(&A, matsize);
    cudaMallocManaged(&B, matsize);
    cudaMallocManaged(&C, matsize);

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            A[N*i+j] = drand48();
            B[N*i+j] = drand48();
        }
    }
    
    double comp_time = 0, comm_time = 0;

    for(int irank=0; irank<mpisize; irank++) { // each rank 
        dim3 grid(N/M, N);
        auto tic = chrono::steady_clock::now();
        //offset = N/mpisize*((mpirank+irank) % mpisize);
        /*matrix multiplication*/
        /***
        for (int i=0; i<N/mpisize; i++)
            for (int j=0; j<N/mpisize; j++)
                for (int k=0; k<N; k++)
                    subC[N*i+j+offset] += subA[N*i+k] * subB[N/mpisize*k+j]; //subC[i, j] = subA[i, k] * subB[k, j]
        ***/  
        matmul<<<grid, M>>>(A, B, C, N); 
        // N/M >= 1 
        cudaDeviceSynchronize();
        auto toc = chrono::steady_clock::now();
        comp_time += chrono::duration<double>(toc - tic).count();
        tic = chrono::steady_clock::now();
        comm_time += chrono::duration<double>(tic - toc).count();
  }

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            for (int k=0; k<N; k++){
                C[N*i+j] -= A[N*i+k] * B[N*k+j];
            }
        }
    }
  
    double err = 0;
    
    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            err += fabs(C[N*i+j]);
        }   
    }
    
    if(mpirank==0) {
        double time = comp_time+comm_time;
        printf("N    : %d\n",N);
        printf("comp : %lf s\n", comp_time);
        printf("comm : %lf s\n", comm_time);
        printf("total: %lf s (%lf GFlops)\n",time,2.*N*N*N/time/1e9);
        printf("error: %lf\n",err/N/N);
    }

    cudaFree(A);
    cudaFree(B);
    cudaFree(C);

    MPI_Finalize();
}
