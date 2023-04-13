import os
import sys
import argparse
import yaml
from datetime import datetime, date

def MC(inputCfg,runmode):
    
    fIn  = open(inputCfg["MC"]["run_list_file"], "r")
    dir=("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s/%s" % (inputCfg["MC"]["task_version"],runmode,inputCfg["MC"]["macro_to_run"]))
    if not os.path.isdir("%s" % (dir)):
        os.system("mkdir -p %s" % (dir))
    
    for run in fIn:
        print(fr"root -b -q %s\(\"%s\",%i,\"%s\",kFALSE,\"%s\")" %(inputCfg["MC"]["macro_to_run"],runmode,int(run),inputCfg["MC"]["task_version"],inputCfg["MC"]["prefix_file_name"]))
        os.system(fr"root -b -q %s.C+\(\"%s\",%i,\"%s\",kFALSE,\"%s\"\)" %(inputCfg["MC"]["macro_to_run"],runmode,int(run),inputCfg["MC"]["task_version"],inputCfg["MC"]["prefix_file_name"]))
def merge_files(inputCfg,runmode):
    '''
    function for merging files
    '''
    dir=("/home/michele_pennisi/cernbox/output_HF_dimuons/mc_analysis_output/%s/%s/%s" % (inputCfg["MC"]["task_version"],runmode,inputCfg["MC"]["macro_to_run"]))
    mergeCommand = "hadd -f -k %s/%s_%s_merged.root " % (dir,runmode,inputCfg["merging_MC"]["input_file"])
    
    fIn  = open(inputCfg["MC"]["run_list_file"], "r")
    for run in fIn:
        partialCommand = "  %s/%s_%s_%i.root" % (dir,runmode,inputCfg["merging_MC"]["input_file"],int(run))
        print(partialCommand,end="\n",flush=True)
        mergeCommand= mergeCommand+partialCommand
        
    print(mergeCommand,end="\n",flush=True)
    print('Proceed with the merging? (true / false)')    
    merge_exec = input()
    if merge_exec == 'true' :
        os.system(mergeCommand)
        


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--mc", help="merge output files", action="store_true")
    parser.add_argument("--merge", help="merge output files", action="store_true")
    parser.add_argument("--HF", help="merge output files", action="store_true")
    parser.add_argument("--MB", help="merge output files", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')
    if args.merge:
        if args.HF:
            merge_files(inputCfg,"HF")
        if args.MB:
            merge_files(inputCfg,"MB")
    if args.mc:
        if args.HF:
            MC(inputCfg,"HF")
        if args.MB:
            MC(inputCfg,"MB")

main()
