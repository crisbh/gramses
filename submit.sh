#!/bin/sh

#SBATCH -n 12
#SBATCH -t 72:00:00
#SBATCH -J gramses
#SBATCH -o ./logs/%j.log
#SBATCH -e ./logs/%j.err
#SBATCH -p cosma
#SBATCH -A durham
#SBATCH --exclusive

# prepare
module purge
module load intel_comp
module load intel_mpi
module load hdfview
unset I_MPI_HYDRA_BOOTSTRAP

# compile
#make clean && make

# run
mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_test.nml

# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode
exit
