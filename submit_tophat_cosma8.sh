#!/bin/sh

#SBATCH -n 128
#SBATCH -t 02:00:00
#SBATCH -J tophat
#SBATCH -o ./logs/%j.log
#SBATCH -e ./logs/%j.err
#SBATCH -p cosma8
#SBATCH -A dp203 
#SBATCH --exclusive

# prepare
module purge
module load intel_comp/2018-update2
module load intel_mpi/2018
# module load hdfview
unset I_MPI_HYDRA_BOOTSTRAP

# compile
cd bin && make clean && make && cd ..

# run
# smooth density field with parts at cell boundaries
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM64_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B256_PM512_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM64_test_cosma8.nml
mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B256_PM128_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM128_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM256_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM256_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM512_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B1243_PM256_cosma8.nml

# Back up simulation data
echo "Job done. Now backing up data..."
JOB_DIR_NAME="data_job_"$SLURM_JOBID
mkdir $JOB_DIR_NAME
cp -r output_000* $JOB_DIR_NAME
cp logs/$SLURM_JOBID.* $JOB_DIR_NAME


# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode
exit
