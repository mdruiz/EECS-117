#!/bin/bash
#$ -N Mandelbrot
#$ -q eecs117,pub64,pub8i
#$ -pe mpi 64
#$ -R y

# Grid Engine Notes:
# -----------------
# 1) Use "-R y" to request job reservation otherwise single 1-core jobs
#    may prevent this multicore MPI job from running.   This is called
#    job starvation.

# Module load boost
module load boost/1.57.0

# Module load OpenMPI
module load openmpi-1.8.3/gcc-4.9.2

# Run the program 
mpirun -np 1  ./mandelbrot_serial 5000 5000

mpirun -np 64  ./mandelbrot_susie 5000 5000

mpirun -np 64  ./mandelbrot_joe 5000 5000

mpirun -np 64  ./mandelbrot_ms 5000 5000



mpirun -np 1  ./mandelbrot_serial 10000 10000

mpirun -np 64  ./mandelbrot_susie 10000 10000

mpirun -np 64  ./mandelbrot_joe 10000 10000

mpirun -np 64  ./mandelbrot_ms 10000 10000