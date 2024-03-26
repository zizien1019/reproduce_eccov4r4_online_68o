#!/bin/bash
# The interpreter used to execute the script

#SBATCH --job-name=MiPla
#SBATCH --mail-user=zizien@umich.edu
#SBATCH --mail-type=BEGIN,END,FAIL
#SBATCH --nodes=3
#SBATCH --ntasks-per-node=32
#SBATCH --cpus-per-task=1
#SBATCH --mem-per-cpu=200m
#SBATCH --time=03:00:00
#SBATCH --account=yulinpan0
#SBATCH --partition=standard

export LANG=en_US.utf8
export LC_ALL=en_US.utf8

# # # Tests Reminder
# ** in Quasi-2D setting
# 1.  I.C.: surface constant     Source: zero
# 1.5 I.C.: gdp                  Source: zero
# 2.  I.C.: zero                 Source: GDP
# 3.  I.C.: zero                 Source: GDP / 100

# # # Step 1 : build
# module purge
# module load intel impi
# cd /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/build
# make CLEAN
# # # for 2d case
# ../../MITgcm/tools/genmake2 -mods=../*68o/zero_w -rd=../../MITgcm -optfile=../*68o/linux_amd64_ifort+impi -mpi
# # # for 3d case
# # ../../MITgcm/tools/genmake2 -mods=../*68o/wrise_rhop -rd=../../MITgcm -optfile=../*68o/linux_amd64_ifort+impi -mpi
# make -j96 depen
# make -j96 all




# # # Step 2 : run
module purge
module load intel
module load impi
module load hdf5
module load netcdf-fortran
unset I_MPI_PMI_LIBRARY
export I_MPI_JOB_RESPECT_PROCESS_PLACEMENT=0

cd /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/run
# rm -rf /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/run/*

# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/input_init/error_weight/data_error/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/input_init/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/input_forcing/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/other/flux-forced/state_weekly/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/other/flux-forced/forcing/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/other/flux-forced/forcing_weekly/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/other/flux-forced/mask/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/forcing/other/flux-forced/xx/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/gcmfaces_climatologies/*.bin .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/ic_files/* .

# most important :: specify input files
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/reproduce_eccov4r4_online_68o/frag_input/* .
ln -s /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/reproduce_eccov4r4_online_68o/zerow_input/* .
# ln -s /scratch/yulinpan_root/yulinpan98/zizien/ecco_v4r4/reproduce_eccov4r4_online_68o/wrise_rhop_input/* .

cp -p ../build/mitgcmuv .
mpiexec -np 96 ./mitgcmuv > a.log

# # Step 3 : Store Data
mkdir ziens_S
# # mv S.0* T.0* U.0* V.0* W.0* ziens_stuvw
mv PTR* ziens_S
mv ziens_S ..