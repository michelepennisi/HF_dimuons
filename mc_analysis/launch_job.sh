#!/bin/bash 

file=$1
readarray -t Arr < ${file}
echo ${#Arr[@]}

mode=${Arr[1]}
task_version=${Arr[2]}
mc_type=${Arr[3]}
grid_work_dir=${Arr[4]}

working_dir=/alidata/mpennisi/HF_dimuons/mc_analysis/analysis_grid/grid_sim/${Arr[0]}

echo mkdir -p ${working_dir} 
echo cd ${working_dir}
echo rm AliAnalysisTaskDimuon_HighMass.* AddTaskDimuon_HighMass.C
echo cp ../AliAnalysisTaskDimuon_HighMass.* ../AddTaskDimuon_HighMass.C .

for (( i=5; i<${#Arr[@]}; i++ )); 

do 
echo aliroot -b -q ../ReadMCDimuon_HighMass.C+\(\"${mode}\",${Arr[i]},\"${task_version}\",\"${mc_type}\",\"${grid_work_dir}\"\)

done