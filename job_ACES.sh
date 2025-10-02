#!/bin/bash
#SBATCH --job-name=MP+BFLG     # Set the job name
#SBATCH --time=3-00:00:00                 # Set the wall clock limit
#SBATCH --nodes=2                           
#SBATCH --ntasks=192
#SBATCH --ntasks-per-node=96
#SBATCH --mem=200G
#SBATCH --output=%x.%j                  # Redirect stdout/stderr to a file
#SBATCH --partition=cpu                 # Specify partition to submit job to
#SBATCH --mail-user=zizien@umich.edu
#SBATCH --mail-type=ALL

# Exit immediately if a command fails, if an unset variable is used,
# and propagate errors in pipelines (safer/stricter mode).
set -euo pipefail








export LANG=en_US.utf8
export LC_ALL=en_US.utf8

RUNNER_HOME='/scratch/user/u.zt173227/'
SIMULA_HOST='/scratch/user/u.zt173227/BFLG/'
ECCOv4r5_DD='/scratch/group/p.phy250235.000/alan/Release5/'

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 1 : build  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%

cd $SIMULA_HOST

rm -rf $SIMULA_HOST/code/*
cp $RUNNER_HOME/git/reproduce_eccov4r4_online_68o/ecco_v4r5/code/*    $SIMULA_HOST/code/
cp $RUNNER_HOME/git/reproduce_eccov4r4_online_68o/i_from_ecco/code/*  $SIMULA_HOST/code


cd $SIMULA_HOST/build

module purge
module load intel-compilers/2021.4.0
module load impi/2021.4.0

make CLEAN
$RUNNER_HOME/MITgcm/tools/genmake2 -mods=$SIMULA_HOST/code -rd=$RUNNER_HOME/MITgcm/ -optfile=$SIMULA_HOST/linux_amd64_ifort+impi -mpi

make depend
make all

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 2 : run  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# arrange files in run/


cd $SIMULA_HOST/run

rm -rf $SIMULA_HOST/run/*

ln -s $RUNNER_HOME/gcmfaces_climatologies/* .
ln -s $ECCOv4r5_DD/freshwater_runoff/* .
ln -s $ECCOv4r5_DD/TBADJ .

rm -rf $SIMULA_HOST/input/*
cp $RUNNER_HOME/git/reproduce_eccov4r4_online_68o/ecco_v4r5/input/*    $SIMULA_HOST/input/
cp $RUNNER_HOME/git/reproduce_eccov4r4_online_68o/i_from_ecco/input/*  $SIMULA_HOST/input
ln -s $SIMULA_HOST/input/* .
ln -s $ECCOv4r5_DD/input_bin/* .

cp -p $SIMULA_HOST/build/mitgcmuv .












# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# set mpi environment and run

module purge
module load intel-compilers/2021.4.0
module load impi/2021.4.0

unset I_MPI_PMI_LIBRARY                           # let Slurm provide PMI
unset I_MPI_PMI                                   # ensure not set to pmix
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0      # avoid binding warnings
export I_MPI_FABRICS=shm:ofi                      # shared memory + OFI
export I_MPI_OFI_PROVIDER=verbs                   # InfiniBand via OFI verbs
export FI_PROVIDER=verbs                          # ensure libfabric uses verbs
# Optional shims if connection issues occur:
# export FI_OFI_RXM_USE_SRX=1                     # RXM shim for scalable EPs
# export FI_VERBS_IFACE=ib0                       # pin to specific IB iface if needed
export I_MPI_DEBUG=5                              # increase IMPI debug verbosity

echo "Nodes=$SLURM_JOB_NUM_NODES  Tasks=$SLURM_NTASKS  NPN=$SLURM_NTASKS_PER_NODE"
echo "MPI launcher: srun --mpi=pmi2"
echo "I_MPI_FABRICS=$I_MPI_FABRICS"
echo "I_MPI_OFI_PROVIDER=$I_MPI_OFI_PROVIDER"
echo "FI_PROVIDER=${FI_PROVIDER:-}"
cd $SIMULA_HOST/run
srun --mpi=pmi2 -n ${SLURM_NTASKS} $SIMULA_HOST/run/mitgcmuv > $SIMULA_HOST/run/a.log

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
########################################%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 3 : post-process  %%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

# cd $SIMULA_HOST
# mkdir output
# mv $SIMULA_HOST/run/PTRACER* $SIMULA_HOST/output/
# cp $SIMULA_HOST/input/data.ptracers $SIMULA_HOST/output/
# cp $SIMULA_HOST/input/data.gchem $SIMULA_HOST/output/
# cp $SIMULA_HOST/code/gchem_calc_tendency.F $SIMULA_HOST/output/

# mv $SIMULA_HOST/run/pickup* $SIMULA_HOST/input/
