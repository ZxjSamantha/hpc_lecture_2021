module load intel-mpi
module load cuda openmpi
module load gcc/8.3.0
module load cuda
nvcc 01_cuda.cu -lmpi
./a.out
