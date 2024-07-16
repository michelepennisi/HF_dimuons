#!/bin/bash

# Execute the actual job script

hostname

#dir_in=/eos/user/m/mpennisi/Data/output_grid/$1/Version_5_AliAOD_skimmed_fwd_fullstat

#echo $dir_in

#ls $dir_in

#ls /eos/user/m/mpennisi/Data/output_grid

#dir_out=/eos/user/m/mpennisi/Data/output_grid/$1/Version_5_AliAOD_skimmed_fwd_fullstat_mod

echo python /afs/cern.ch/work/m/mpennisi/private/powheg_std/run_powheg_sim.py "$@"

python3 /afs/cern.ch/work/m/mpennisi/private/powheg_std/run_powheg_sim.py "$@"


pwd