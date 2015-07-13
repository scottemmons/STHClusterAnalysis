## Author: Scott Emmons (scott@scottemmons.com)
## Purpose: To automate the command line calls necessary to execute the full clustering analysis script workflow.
## Date: January 2, 2015

## Depends on compiled LFR benchmark graph generation code and compiled Lancichinetti clustering repository code.

import argparse
import os
import shutil
import errno
import time 
import math
import random
import glob
import subprocess

supported_modules = ['generate', 'cluster', 'measure', 'visualize']
bench_directory_lock = 'bench_directory.lock'
lancichinetti_program_lock = 'lancichinetti_program.lock'

def handleArgs():
    """Handles the command-line input arguments, placing them in the global Namespace variable 'args'."""

    parser = argparse.ArgumentParser(description="Automates the command line calls necessary to execute the full clustering analysis script workflow")
    parser.add_argument("--mu", type=float, required=True, help="the mixing parameter. Fortunato suggests we test values of 0.40, 0.50, and 0.60", dest="mu")
    parser.add_argument("-t", "--times", default=1, type=int, help="the number of clusterings to run on each graph", dest="numtrials")
    parser.add_argument("-n", "--nlist", nargs="+", required=True, help="the values of N, or number of nodes in the graph, for which to execute the workflow", dest="n_list")
    parser.add_argument("-s", "--start", default = 1, type=int, help="the file number at which to start, inclusive", dest="start")
    parser.add_argument("-e", "--end", default = 10, type=int, help="the file number at which to end, inclusive", dest="end")
    parser.add_argument("-m", "--modules", nargs="+", default=supported_modules, type=str.lower, choices=supported_modules, help="the names of the modules to run. defaults to all of them", dest="modules")
    parser.add_argument("-b", "--benchmark", default=os.getcwd()+"/binary_networks/", help="the path to the installed LFR generation software, required for functionality of the generate module. Defaults to the current working directory + '/binary_networks/'", dest="bench_directory_stem")
    parser.add_argument("--lp", default=os.getcwd()+"/clustering_programs_5_2/", help="the path at which the Lancichinetti clustering program is installed, required for full functionality of the measure module. Defaults to the current working directory + '/clustering_programs_5_2/'", dest="lancichinetti_program_path")

    global args
    args = parser.parse_args()

def touch(filename, times=None):
    """Python implementation of 'touch' command in Ubuntu terminal. Not completely race-free. Credits to community wiki on SO: def touch(fname, times=None): http://stackoverflow.com/questions/1158076/implement-touch-using-python"""
    with open(filename, 'a'):
        os.utime(filename, times)

def combineDataFiles(name, extension, directory):
    """Combines the data files matching directory + '/' + name + '_s*_e*' + extension in directory.
    The arguments name, extension, and directory should be strings."""

    def deleteFileIfNeeded(file):
        try:
            os.remove(file)
        except OSError as error:
            if error.errno != errno.ENOENT:
                raise

    def appendFileToFile(append_file_name, to_file_name, write_header = True):
        """Write the content in append_file_name into to_file_name. Write header determines whether or not the first line in append_file_name is written to to_file_name."""

        to_f = open(to_file_name, 'a')
        append_f = open(append_file_name, 'r')

        if not write_header:
            append_f.readline()
        for line in append_f:
            to_f.write(line)

        to_f.close()
        append_f.close()

    filtered_files = glob.glob(directory + '/' + name + '_s*_e*' + extension)
    assert len(filtered_files) > 0

    combined_file = directory + '/' + name + extension
    deleteFileIfNeeded(combined_file)
    written_header = False
    for data_file in filtered_files:
        appendFileToFile(data_file, combined_file, not written_header)
        written_header = True

def incrementBenchSeed(bench_directory_stem, times = 1):
    """Increments the value of the random seed in Lancichinetti's bencmark generation program."""

    seed_file = bench_directory_stem + 'time_seed.dat'
    
    with open(seed_file, 'r') as f:
        value = f.readline().split()[0]

    with open(seed_file, 'wb') as f:
        f.write(str(int(value) + 1 * times) + '\n')

handleArgs()

#To scatter parallel runs executed simultaneously
##time.sleep(random.randrange(500, 18000) / 100.0)

##args.n_list = ['10', '100', '1000', '10000', '100000']

for n in args.n_list:

    working_directory = 'generated_benches/n_' + n + '/'
    new_bench_directory_stem = working_directory + 'binary_networks_s_' + str(args.start) + '_e_' + str(args.end) + '/'
    new_lancichinetti_program_path = working_directory + 'clustering_programs_5_2_s_' + str(args.start) + '_e_' + str(args.end) + '/'

    if 'generate' in args.modules:

        # python generate.py -n 1000 --mu 0.40 -s 1 -e 10 -b binary_networks/ -o generated_benches/n_1000/
        command_list = []
        command_list.append('python')
        command_list.append('generate.py')
        command_list.append('-n')
        command_list.append(n)
        command_list.append('--maxk')
        command_list.append(str(int(0.10 * float(n))))
        command_list.append('--mu')
        command_list.append(str(args.mu))
        command_list.append('--minc')
        command_list.append('50')
        command_list.append('--maxc')
        command_list.append(str(int(0.10 * float(n))))
        command_list.append('-s')
        command_list.append(str(args.start))
        command_list.append('-e')
        command_list.append(str(args.end))
        command_list.append('-b')
        command_list.append(new_bench_directory_stem)
        command_list.append('-o')
        command_list.append(working_directory)

        #To be run in parallel
	bench_copy_start = time.clock()
        shutil.copytree(args.bench_directory_stem, new_bench_directory_stem)
	bench_copy_end = time.clock()
	print '\nTime to copy bench_directory: ' + str(bench_copy_end - bench_copy_start) + '\n'
        increment_by = args.start + (int(math.ceil(math.log(int(n), 10))) - 3) * 100
        incrementBenchSeed(new_bench_directory_stem, times = increment_by)
        print '\nIncremented bench seed by: ' + str(increment_by) + '\n'
        touch('timeseed_inc_token_n_' + n + '_s_' + str(args.start) + '_e_' + str(args.end))

        subprocess.call(command_list)
        bench_remove_start = time.clock()
        shutil.rmtree(new_bench_directory_stem)
        bench_remove_end = time.clock()
        print '\nTime to remove bench_directory: ' + str(bench_remove_end - bench_remove_start) + '\n'

        print '\n' + ('*' * (62 + len(n)))
        print '* Completed clustering analysis "generate" workflow for N = ' + n + ' *'
        print '*' * (62 + len(n)) + '\n'

    if 'cluster' in args.modules:
        
        # python cluster.py -m blondel infomap label_propagation slm -t 1 -u --gpre generated_benches/n_1000/network_v --gsuf .dat -s 1 -e 10 --lp clustering_programs_5_2/ -o generated_benches/n_1000/
        # python cluster.py -m blondel infomap label_propagation slm -t 1 -u --gpre real_world_graphs/flickrEdges_renum_v --gsuf .txt -s 1 -e 1 --lp clustering_programs_5_2/ -o real_world_graphs/
        command_list = []
        command_list.append('python')
        command_list.append('cluster.py')
        command_list.append('-m')
        command_list.append('blondel')
        command_list.append('infomap')
        command_list.append('label_propagation')
        command_list.append('slm')
        command_list.append('-t')
        command_list.append(str(args.numtrials))
        command_list.append('-u')
        command_list.append('--gpre')
        command_list.append(working_directory + 'network_v')
        command_list.append('--gsuf')
        command_list.append('.dat')
        command_list.append('-s')
        command_list.append(str(args.start))
        command_list.append('-e')
        command_list.append(str(args.end))
        command_list.append('--lp')
        command_list.append(new_lancichinetti_program_path)
        command_list.append('-o')
        command_list.append(working_directory)

        #To be run in parallel
	lancichinetti_copy_start = time.clock()
        shutil.copytree(args.lancichinetti_program_path, new_lancichinetti_program_path)
	lancichinetti_copy_end = time.clock()
	print '\nTime to copy lancichinetti_program: ' + str(lancichinetti_copy_end - lancichinetti_copy_start) + '\n'

        subprocess.call(command_list)
        lancichinetti_remove_start = time.clock()
        shutil.rmtree(new_lancichinetti_program_path)
        lancichinetti_remove_end = time.clock()
        print '\nTime to remove lancichinetti_program: ' + str(lancichinetti_remove_end - lancichinetti_remove_start) + '\n'

        print '\n' + ('*' * (61 + len(n)))
        print '* Completed clustering analysis "cluster" workflow for N = ' + n + ' *'
        print '*' * (61 + len(n)) + '\n'

    if 'measure' in args.modules:
        
        # python measure.py --gmap gmap/ --lnmi mutual3/ --gpre generated_benches/n_1000/network_v --gsuf .dat --gsep $'\t' -u --srun --spre generated_benches/n_1000/community_v --ssuf .dat --ssep $'\t' -s 1 -e 10 --cnames blondel infomap label_propagation slm --cpre generated_benches/n_1000/blondel_clustering_v generated_benches/n_1000/infomap_clustering_v generated_benches/n_1000/label_propagation_clustering_v generated_benches/n_1000/slm_clustering_v --csuf .dat --csep $'\t' --cnum 10 -o generated_benches/n_1000/
        command_list = []
        command_list.append('python')
        command_list.append('measure.py')
        command_list.append('--gmap')
        command_list.append('gmap/')
        command_list.append('--lnmi')
        command_list.append('mutual3/')
        command_list.append('--gpre')
        command_list.append(working_directory + 'network_v')
        command_list.append('--gsuf')
        command_list.append('.dat')
        command_list.append('--gsep')
        command_list.append('\t')
        command_list.append('-u')
        command_list.append('--srun')
        command_list.append('--spre')
        command_list.append(working_directory + 'community_v')
        command_list.append('--ssuf')
        command_list.append('.dat')
        command_list.append('--ssep')
        command_list.append('\t')
        command_list.append('-s')
        command_list.append(str(args.start))
        command_list.append('-e')
        command_list.append(str(args.end))
        command_list.append('--cnames')
        command_list.append('blondel')
        command_list.append('infomap')
        command_list.append('label_propagation')
        command_list.append('slm')
        command_list.append('--cpre')
        command_list.append(working_directory + 'blondel_clustering_v')
        command_list.append(working_directory + 'infomap_clustering_v')
        command_list.append(working_directory + 'label_propagation_clustering_v')
        command_list.append(working_directory + 'slm_clustering_v')
        command_list.append('--csuf')
        command_list.append('.dat')
        command_list.append('--csep')
        command_list.append('\t')
        command_list.append('--cnum')
        command_list.append(str(args.numtrials))
        command_list.append('-o')
        command_list.append(working_directory)
        subprocess.call(command_list)

        print '\n' + ('*' * (61 + len(n)))
        print '* Completed clustering analysis "measure" workflow for N = ' + n + ' *'
        print '*' * (61 + len(n)) + '\n'

if 'visualize' in args.modules:

    # Combine output raw data from previous step
    for n in args.n_list:
        combineDataFiles('raw_data', '.csv', 'generated_benches/n_' + n)

    # Combined output timing data from previous step
    for n in args.n_list:
        combineDataFiles('timing_data', '.csv', 'generated_benches/n_' + n)

    # python visualize.py -f generated_benches/n_100/raw_data.csv generated_benches/n_1000/raw_data.csv generated_benches/n_10000/raw_data.csv --names n_100 n_1000 n_10000 --nodes 100 1000 10000
    command_list = []
    command_list.append('python')
    command_list.append('visualize.py')
    command_list.append('-f')
    for n in args.n_list:
        command_list.append('generated_benches/n_' + n + '/raw_data.csv')
    command_list.append('--names')
    for n in args.n_list:
        command_list.append('n_' + n)
    command_list.append('--nodes')
    for n in args.n_list:
        command_list.append(n)
    subprocess.call(command_list)

    print '\n' + ('*' * (63 + len(n)))
    print '* Completed clustering analysis "visualize" workflow for N = ' + n + ' *'
    print '*' * (63 + len(n)) + '\n'
