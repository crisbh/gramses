#!/bin/sh

#SBATCH -n 256
#SBATCH -t 02:00:00
#SBATCH -J tophat
#SBATCH -o ./logs/%j.log
#SBATCH -e ./logs/%j.err
#SBATCH -p cosma7
#SBATCH -A dp004 
#SBATCH --exclusive

# prepare
module purge
module load intel_comp/2018-update2
module load intel_mpi/2018
# module load hdfview
unset I_MPI_HYDRA_BOOTSTRAP

# compile
#make clean && make

# run
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B512_PM256.nml
mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM256.nml

# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode
exit