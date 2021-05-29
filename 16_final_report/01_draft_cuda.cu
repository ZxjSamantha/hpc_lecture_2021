#include <cstdio> 
#include <cstdlib>
#include <cmath>
#include <vector>
#include <chrono> 
using namespace std;

__global__ void subA(float *out, float *a, float *b, int N, int offset)
{
    for (int i=0; i < N; i++)
        for (int j=0; j < N; j++)
            out[N*i+j] = A[N*(i+offset)+j]    
}

__global__ void subB(float *out, float *a, float *b, int N, int offset, int )
{
    for (int i=0; i < N; i++)
        for (int j=0; j < N; j++)
            out[N/size*i+j] = B[N*i+j+offset]
}

int main(int argc, char** argv){
    int size, rank; 
    MPI_Init(&argc, &argv); 
    MPI_Comm_size(MPI_COMM_WORLD, &size); 
    MPI_Comm_size(MPI_COMM_WORLD, &rank); 

    const int N = 256; 

    vector<float> A(N*N);
    vector<float> B(N*N);
    vector<float> C(N*N, 0);
    vector<float> subA(N*N/size); 
    vector<float> subB(N*N/size);
    vector<float> subC(N*N/size, 0);

    for (int i=0; i<N; i++){
        for (int j=0; j<N; j++){
            A[N*i+j] = drand48();
            B[N*i+j] = drand48();
        }
    }
    int offset = N/size*rank; 

    for (int i=0; i<N/size; i++)
        for (int j=0; j<N; j++)
            subA[N*i+j] = A[N*(i+offset)+j]; 
    for (int i=0; i<N/size; i++)
        for (int j=0; j<N; j++)
            subB[N*i+j] = B[N*(i+offset)+j]; 
}


