import sys
import subprocess
import time
import shutil
import os

def main(args):
    print("entro nel main")
    os.system("hostname")
    # os.system("/cvmfs/alice.cern.ch/bin/alienv enter VO_ALICE@AliPhysics::v5-09-56t-01_O2-1,VO_ALICE@AliDPG::prod-202307-02-1,VO_ALICE@POWHEG::r3964-alice2-5")
    output_file = "output.txt"
    # Example usage:
    file_path = 'powheg.input'  # Replace with your file path
    
    replacements = create_dictionary(args)
    replace_in_file(file_path, replacements)

    #cluster_id = os.getenv('CONDOR_PROCNO')
    #cluster_id = "001"
    
    command = "pwhg_main_Z"

    # command = "ls"
    
    try:
        process = subprocess.Popen(command, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
        stdout_bytes, stderr_bytes = process.communicate()
        
        # Decode bytes to strings
        stdout_str = stdout_bytes.decode('utf-8')
        stderr_str = stderr_bytes.decode('utf-8')
        
        # Print output
        print("STDOUT:\n", stdout_str)
        print("STDERR:\n", stderr_str)
        # Write command output to output_file if needed
        with open(output_file, "w") as f:
            f.write(stdout_str)
            f.write("Arguments: %s\n" % args)

    except subprocess.CalledProcessError as e:
        print("Error: %s\n" % {e.stderr})
        print("Return Code: %s\n" % {e.returncode})
        
    destination_dir = "/afs/cern.ch/work/m/mpennisi/private/powheg_std/output/mass_range_%s_%s/%s" % (args[1], args[2], args[3])

    if not os.path.exists(destination_dir):
        os.makedirs(destination_dir)

    os.system('cp * %s' % (destination_dir))

def replace_in_file(file_path, replacements):
    # Read the entire file
    with open("base_powheg_DY.input", 'r') as file:
        file_content = file.read()

    # Perform replacements
    for old_string, new_string in replacements.items():
        file_content = file_content.replace(old_string, new_string)

    # Write back to the file
    with open(file_path, 'w') as file:
        file.write(file_content)

def create_dictionary(args):
    # Initialize an empty dictionary
    arguments_dict = {
        'numevts 12345': 'numevts %s' % args[0],
        'iseed    654321': 'iseed    %d' % (int(time.time())+int(args[3])),
        'mass_low 4': 'mass_low %s' % args[1],
        'mass_high 40': 'mass_high %s' % args[2]
    }
    
    return arguments_dict

if __name__ == "__main__":
    main(sys.argv[1:])
