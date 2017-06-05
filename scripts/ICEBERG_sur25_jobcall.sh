#!/bin/bash
#$ -l h_cpu=70:00:00
#$ -pe openmp 4
#$ -l mem=2G
#$ -l rmem=2G
#$ -M shaun.coutts@gmail.com
#$ -m be
#$ -cwd

module load apps/R/3.3.1
module add mpi/gcc/openmpi/1.8.8
mpirun -np 1 Rscript /home/bo1src/ind_dem_perf/scripts/SPP_25_sur_stan.R > /home/bo1src/ind_dem_perf/output/console/SPP_25_sur_cons.out 2>&1
