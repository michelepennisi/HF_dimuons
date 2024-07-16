#!/bin/bash

# Execute the actual job script

hostname

dir=/afs/cern.ch/work/m/mpennisi/private/powheg_std/output/mass_range_$2_$3/$4

python /afs/cern.ch/work/m/mpennisi/private/powheg_std/LHE_NTuples/LHEConverter.py -i $dir/pwgevents.lhe -o $dir/converted_lhe.root
