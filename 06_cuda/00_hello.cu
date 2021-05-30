#include <cstdio>

__global__ void print(void) { // Call GPU on CPU
  printf("Hello GPU\n");
}

int main() {
  printf("Hello CPU\n");
  //print<<<1,1>>>();//Change the number 
  //foo<<<blocks, threads_per_block>>>(); 
  print<<<1,4>>>();//2048 wont work; 1024 is OK 
  cudaDeviceSynchronize();
}
