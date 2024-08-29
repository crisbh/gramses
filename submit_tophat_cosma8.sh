#!/bin/sh

#SBATCH -n 1024
#SBATCH -t 32:00:00
#SBATCH -J 9e-2AMR
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

# Define input param file
PARAM_FILE="./namelist/gr_tophat_B256_PM512_deltai_9e-2_AMR_cosma8.nml"
#PARAM_FILE="./namelist/gr_tophat_B256_PM512_deltai_9e-2_cosma8.nml"
#PARAM_FILE="./namelist/gr_tophat_B256_PM512_cosma8.nml"
#PARAM_FILE="./namelist/gr_tophat_B25_6_PM512_cosma8.nml"

# run
# smooth density field with parts at cell boundaries
mpirun -np $SLURM_NTASKS ./bin/ramses3d $PARAM_FILE
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM256_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM512_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B256_PM512_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B256_PM128_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM128_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM256_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_PM512_cosma8.nml
#mpirun -np $SLURM_NTASKS ./bin/ramses3d ./namelist/gr_tophat_B1243_PM256_cosma8.nml

# Post processing
## Readgrav data
#bash ./run_batch_readgrav.sh $SLURM_NTASKS 9
#
## Back up simulation data
#echo "Job done. Now backing up data..."
#JOB_DIR_NAME=snapshots/"data_job_"$SLURM_JOBID
#mkdir $JOB_DIR_NAME
#mv output_000* $JOB_DIR_NAME
#cp logs/$SLURM_JOBID.* $JOB_DIR_NAME


# report
echo ""
sacct -j $SLURM_JOBID --format=JobID,JobName,Partition,AllocCPUS,Elapsed,ExitCode
exit
