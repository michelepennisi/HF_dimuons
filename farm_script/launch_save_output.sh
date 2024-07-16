#!/bin/bash
# filename=$1
# while read -r period run_number  gen ; do
#   echo $period
#   
#   echo ${dir_in}
#   echo $run_number
#   echo ${gen}
  
  
# done < $filename
hostname

echo "Entering the main script changed 6"

# # Check if CVMFS is mounted and source the environment
# if [ -d /cvmfs/alice.cern.ch ]; then
#     echo "CVMFS is mounted. Attempting to source the environment."
#     source /cvmfs/alice.cern.ch/etc/login.sh
#     echo /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@AliPhysics::v5-09-56t-01_O2-1,VO_ALICE@AliDPG::prod-202307-02-1,VO_ALICE@POWHEG::r3964-alice2-5
#     /cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@AliPhysics::v5-09-56t-01_O2-1,VO_ALICE@AliDPG::prod-202307-02-1,VO_ALICE@POWHEG::r3964-alice2-5
#     echo alienv list
#     alienv list
# else
#     echo "CVMFS is not mounted. Exiting."
#     exit 1
# fi

dir_in=/alidata/mpennisi/output_grid/$1/Version_5_AliAOD_skimmed_fwd_fullstat/$2/output/000 

dir_out=/alidata/mpennisi/output_grid/$1/Version_5_AliAOD_skimmed_fwd_fullstat_mod/$2/output/000 

mkdir -p $dir_out

ls $dir_in

echo root -q /home/mpennisi/HF_dimuons/mc_analysis/read_output/save_mc_output.C\(\"${1}\",\"${dir_in}\",\"${dir_out}\",$2,\"${3}\"\)


root -q /home/mpennisi/HF_dimuons/mc_analysis/read_output/save_mc_output.C\(\"${1}\",\"${dir_in}\",\"${dir_out}\",$2,\"${3}\"\)

#  