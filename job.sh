#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=MiPla
#SBATCH --mail-user=zizien@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=96
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200m
#SBATCH --time=10:00:00
#SBATCH --account=yulinpan0
#SBATCH --partition=standard

export LANG=en_US.utf8
export LC_ALL=en_US.utf8



# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%  Step 1 : build  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%

cd /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/build

module purge
module load intel impi
make CLEAN

# %%%%%%%%%%%%%%%%%%%%%
# # # for 2d case

../../MITgcm/tools/genmake2 -mods=../reproduce_eccov4r4_online_68o/zero_w -rd=../../MITgcm -optfile=../reproduce_eccov4r4_online_68o/linux_amd64_ifort+impi -mpi

# # # for 3d case

# ../../MITgcm/tools/genmake2 -mods=../*68o/wrise_rhop -rd=../../MITgcm -optfile=../*68o/linux_amd64_ifort+impi -mpi

make -j96 depen
make -j96 all




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

cd /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/run

rm -rf ../run/*

ln -s ../../forcing/input_init/error_weight/data_error/* .
ln -s ../../forcing/input_init/error_weight/data_error/* .
ln -s ../../forcing/input_init/* .
ln -s ../../forcing/input_forcing/* .
ln -s ../../forcing/other/flux-forced/state_weekly/* .
ln -s ../../forcing/other/flux-forced/forcing/* .
ln -s ../../forcing/other/flux-forced/forcing_weekly/* .
ln -s ../../forcing/other/flux-forced/mask/* .
ln -s ../../forcing/other/flux-forced/xx/* .
ln -s ../../gcmfaces_climatologies/*.bin .
ln -s ../../reproduce_eccov4r4_online_68o/ic_files/* .

# %%%%%%%%%%%%%%%%%%%%%
# # # most important :: specify input files :: choose one

ln -s ../reproduce_eccov4r4_online_68o/zerow_input/* .
# ln -s ../reproduce_eccov4r4_online_68o/wrise_rhop_input/* .

cp -p ../build/mitgcmuv .
mpiexec -np 96 ./mitgcmuv > a.log




# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%  Step 3 : Store Data  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
# %%%%%%%%%%%%%%%%%%%%%
# # # S.0* T.0* U.0* V.0* W.0* PTR* ziens_stuvw

mkdir ziens_S
mv PTR* ziens_S
mv ziens_S ..