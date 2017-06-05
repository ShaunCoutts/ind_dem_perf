#!/bin/bash
#$ -l h_rt=48:00:00
#$ -pe openmp 4
#$ -l mem=5G
#$ -l rmem=5G
#$ -M shaun.coutts@gmail.com
#$ -m be
#$ -cwd

module load apps/R/3.3.1
module add mpi/gcc/openmpi/1.8.8
mpirun -np 1 Rscript /home/bo1src/ind_dem_perf/scripts/SPP_NK_sur_stan.R > /home/bo1src/ind_dem_perf/output/console/SPP_NK_sur_cons.out 2>&1
