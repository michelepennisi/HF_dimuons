import os
import sys
import argparse
import yaml

def run_on_grid(inputCfg, mode, MC_type):
    '''
    function for run jobs on alien grid
    '''
    if not os.path.isdir(inputCfg["output"]["output_dir"]) :
        print("the directory does not exist, creating %s" % (inputCfg["output"]["output_dir"]))
        os.system("mkdir -p %s" % (inputCfg["output"]["output_dir"]))

    if mode == "test" :
        fIn  = open(inputCfg["input"]["run_list_file"], "r")
        fOut = open("%s/%s_%s_jobs.txt" % (inputCfg["output"]["output_dir"], inputCfg["output"]["submitted_output_file"],MC_type), "w")

    if mode == "full" :
        fIn  = open(inputCfg["input"]["run_list_file"], "r")
        fOut = open("%s/%s_%s_jobs.txt" % (inputCfg["output"]["output_dir"], inputCfg["output"]["submitted_output_file"],MC_type), "w")

    if mode == "terminate" :
        if not os.path.isfile("%s/%s_%s_jobs.txt" % (inputCfg["output"]["output_dir"], inputCfg["output"]["submitted_output_file"],MC_type)) :
            print('Submitted jobs file does not exist! Do --full before')
            return
        fIn  = open("%s/%s_%s_jobs.txt" % (inputCfg["output"]["output_dir"], inputCfg["output"]["submitted_output_file"],MC_type), "r")
        fOut = open("%s/%s_%s_jobs.txt" % (inputCfg["output"]["output_dir"], inputCfg["output"]["terminated_output_file"],MC_type), "w")

    for run in fIn:
        print("aliroot -b -q %s/%s\(\"%s\",%i,\"%s\"\)" % (inputCfg["input"]["dir_macro"],inputCfg["input"]["macro_to_run"],mode, int(run),MC_type))
        os.chdir("%s" % (inputCfg["input"]["dir_macro"]))
        os.system(fr"aliroot -b -q %s\(\"%s\",%i\,\"%s\"\)" % (inputCfg["input"]["macro_to_run"],mode, int(run),MC_type))
        fOut.write(run)

def data(inputCfg):
    
    fIn  = open(inputCfg["input"]["run_list_file"], "r")

    string=("%s/%s" % (inputCfg["data"]["prefix_dir_fileout"],inputCfg["data"]["dir_fileout"]))
    if not os.path.isdir(string):
        os.system("mkdir %s" % (inputCfg["data"]["prefix_dir_fileout"]))
        os.system("mkdir %s/%s" % (inputCfg["data"]["prefix_dir_fileout"],inputCfg["data"]["dir_fileout"]))

    for run in fIn:
          #print(fr"root -b -q %s/%s\(4,%i,\"%s\"\)" % (inputCfg["data"]["macro_dir"],inputCfg["data"]["macro_to_run"],int(run),inputCfg["data"]["dir_fileout"]))
        os.system(fr"root -b -q %s/%s\(4,%i,\"%s\"\)" % (inputCfg["data"]["macro_dir"],inputCfg["data"]["macro_to_run"], int(run),string))

def MC(inputCfg,runmode):
    
    fIn  = open(inputCfg["input"]["run_list_file"], "r")
    if not os.path.isdir("%s/%s" % (inputCfg["MC"]["output_dir"],runmode)):
        os.system("mkdir %s" % (inputCfg["MC"]["output_dir"]))
        os.system("mkdir %s/%s" % (inputCfg["MC"]["output_dir"],runmode))

    for run in fIn:
        print(fr"root -b -q %s/%s\(\"%s\",%i,\"%s\",\"%s\",\"%s\"\)" % (inputCfg["MC"]["macro_dir"],inputCfg["MC"]["macro_to_run"],runmode,int(run),inputCfg["MC"]["input_dir"],inputCfg["MC"]["output_dir"],inputCfg["MC"]["prefix_filename"]))
        os.system(fr"root -b -q %s/%s\(\"%s\",%i,\"%s\",\"%s\",\"%s\"\)" % (inputCfg["MC"]["macro_dir"],inputCfg["MC"]["macro_to_run"],runmode,int(run),inputCfg["MC"]["input_dir"],inputCfg["MC"]["output_dir"],inputCfg["MC"]["prefix_filename"]))

def copy_from_grid(inputCfg):

    if not os.path.isdir(inputCfg["copy"]["output_dir"]) :
        print("the directory does not exist, creating %s" % (inputCfg["output"]["output_dir"]))
        os.system("mkdir -p %s" % (inputCfg["copy"]["output_dir"]))
    '''
    function for downloading files from alien grid
    '''
    if not os.path.isfile("%s/%s" % (inputCfg["copy"]["task_dir"], inputCfg["copy"]["terminated_output_file"])) :
        print('Terminated jobs file does not exist! Do --terminate before')
        return
    fIn  = open("%s/%s" % (inputCfg["copy"]["task_dir"], inputCfg["copy"]["terminated_output_file"]), "r")
    for run in fIn:
        print("alien_cp alien:%s/%i/000/%s_%i.root file:%s/%s_%i.root" % (inputCfg["copy"]["alien_output_path"], int(run),inputCfg["copy"]["file_name"],int(run),inputCfg["copy"]["output_dir"],inputCfg["copy"]["file_name"],int(run)))
        os.system("alien_cp alien:%s/%i/000/%s_%i.root file:%s/%s_%i.root" % (inputCfg["copy"]["alien_output_path"], int(run),inputCfg["copy"]["file_name"],int(run),inputCfg["copy"]["output_dir"],inputCfg["copy"]["file_name"],int(run))) # enable if you want to run on grid

def merge_files(inputCfg,input_type,MC_type):
    '''
    function for merging files
    '''
    if input_type == "Data" :
        
        mergeCommand = "hadd -f %s/%s.root " % (inputCfg["merging_data"]["input_path"], inputCfg["merging_data"]["output_file_prefix"])
        fIn  = open(inputCfg["input"]["run_list_file"], "r")
        for run in fIn:
            mergeCommand += "%s/%s_%i.root " % (inputCfg["merging_data"]["input_path"], inputCfg["merging_data"]["input_file_prefix"], int(run))
            
        
        print(mergeCommand)
        print('Proceed with the merging? (true / false)')    
        merge_exec = input()
        if merge_exec == 'true' :
            os.system(mergeCommand)

    if input_type == "MC":
        mergeCommand = "hadd -f %s/%s/%s_merged.root " % (inputCfg["merging_MC"]["output_path"],MC_type,inputCfg["merging_MC"]["output_file_prefix"])
        
        fIn  = open(inputCfg["input"]["run_list_file"], "r")
        for run in fIn:
            mergeCommand += "%s/%s/%s_%i.root " % (inputCfg["merging_MC"]["output_path"],MC_type,inputCfg["merging_MC"]["input_file"],int(run))
            
        print(mergeCommand)
        print('Proceed with the merging? (true / false)')    
        merge_exec = input()
        if merge_exec == 'true' :
            os.system(mergeCommand)


def main():
    parser = argparse.ArgumentParser(description='Arguments to pass')
    parser.add_argument('cfgFileName', metavar='text', default='config.yml', help='config file name')
    parser.add_argument("--test", help="submit your task in full mode", action="store_true")
    parser.add_argument("--full", help="submit your task in full mode", action="store_true")
    parser.add_argument("--terminate", help="terminate your task", action="store_true")
    parser.add_argument("--copy", help="download files from alien grid", action="store_true")
    parser.add_argument("--mergedata", help="merge output files", action="store_true")
    parser.add_argument("--mergemc", help="merge output files", action="store_true")
    parser.add_argument("--HF", help="merge output files", action="store_true")
    parser.add_argument("--MB", help="merge output files", action="store_true")
    parser.add_argument("--data", help="merge output files", action="store_true")
    parser.add_argument("--MC", help="merge output files", action="store_true")
    args = parser.parse_args()

    print('Loading task configuration: ...', end='\r')
    with open(args.cfgFileName, 'r') as ymlCfgFile:
        inputCfg = yaml.load(ymlCfgFile, yaml.FullLoader)
    print('Loading task configuration: Done!')
    if args.test:
        if args.HF:
            run_on_grid(inputCfg, "test","HF")
        if args.MB:
            run_on_grid(inputCfg, "test","MB")
    if args.full:
        if args.HF:
            run_on_grid(inputCfg, "full","HF")
        if args.MB:
            run_on_grid(inputCfg, "full","MB")
    if args.terminate:
        if args.HF:
            run_on_grid(inputCfg, "terminate","HF")
        if args.MB:
            run_on_grid(inputCfg, "terminate","MB")
    if args.data:
        data(inputCfg)
    if args.copy:
        copy_from_grid(inputCfg)
    if args.mergedata:
        merge_files(inputCfg,"Data","HF")
    if args.mergemc:
        if args.HF:
            merge_files(inputCfg,"MC","HF")
        if args.MB:
            merge_files(inputCfg,"MC","MB")
    if args.MC:
        if args.HF:
            MC(inputCfg,"HF")
        if args.MB:
            MC(inputCfg,"MB")

main()
