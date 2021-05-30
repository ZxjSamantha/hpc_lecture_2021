module load intel-mpi
module load gcc
mpicxx 00_openmp+simd+mpi.cpp -fopt-info-vec-optimized -march=native -fopenmp -O3
mpirun -np 1 ./a.out