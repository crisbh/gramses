#!/bin/sh

#SBATCH -n 128
#SBATCH -t 00:30:00
#SBATCH -J TestTH
#SBATCH -o ./logs/%j.log
#SBATCH -e ./logs/%j.err
#SBATCH -p cosma8
#SBATCH -A dp203 
#SBATCH --exclusive

# prepare
module purge
module load intel_comp
module load compiler-rt tbb compiler mpi
module load python

# compile
cd bin && make clean && make && cd ..

# run
mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B64_PM64.nml

# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode
exit
