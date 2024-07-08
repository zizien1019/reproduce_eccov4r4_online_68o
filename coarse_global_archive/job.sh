#!/bin/bash
#SBATCH --job-name=MiPla
#SBATCH --mail-user=zizien@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=4
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=7g
#SBATCH --time=01:00:00
#SBATCH --account=yulinpan0
#SBATCH --partition=debug

export LANG=en_US.utf8
export LC_ALL=en_US.utf8

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 1 : build  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%

cd /scratch/yulinpan_root/yulinpan98/zizien/tutorial_global_oce_latlon_next/build

module purge
module load intel
module load impi
module load netcdf-fortran
make CLEAN

# %%%%%%%%%%%%%%%%%%%%% make
../../MITgcm/tools/genmake2 -rd=../../MITgcm -mods=../code -optfile=../../ecco_v4r4/reproduce_eccov4r4_online_68o/linux_amd64_ifort+impi -mpi

make -j4 depend
make -j4 all

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 2 : run  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

module purge
module load intel
module load impi
module load hdf5
module load netcdf-fortran
unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

cd /scratch/yulinpan_root/yulinpan98/zizien/tutorial_global_oce_latlon_next/run

rm -rf /scratch/yulinpan_root/yulinpan98/zizien/tutorial_global_oce_latlon_next/run/*

# %%%%%%%%%%%%%%%%%%%%%
ln -s ../input/* .
ln -s ../../ecco_v4r4/reproduce_eccov4r4_online_68o/ic_files/* .
cp -p ../build/mitgcmuv .
mpirun -np 4 ./mitgcmuv > a.log

# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%  Step 3 : Store Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%
# # S.0* T.0* U.0* V.0* W.0* PTR* ziens_stuvw

mkdir ziens_900
mv PTR* ziens_900
mv ziens_900 /scratch/yulinpan_root/yulinpan98/zizien