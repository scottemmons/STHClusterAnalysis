## Author: Scott Emmons
## Purpose: To run various clustering algorithms over the input network graph files.
## Date: December 31, 2014

## Note to self:
## Flesh out parameter and function descriptions

## Current testing command:
## python cluster.py -m demon -t 1 -u --gpre generated_benches/n_1000/network_v --gsuf .dat -s 1 -e 1 --lp clustering_programs_5_2/ -o generated_benches/n_1000/

import argparse
import os
import errno
import shutil
import subprocess
import random
import csv
import datetime

####################
# Global Variables #
####################

supported_methods = ["blondel", "louvain", "label_propagation", "modularity_optimization", "oslom", "infomap", "hierarchical_infomap", "slm", "demon"]
lancichinetti_num_to_name = {0 : 'oslom', 1 : 'oslom', 2 : 'infomap', 3 : 'infomap', 4 : 'blondel', 5 : 'label_propagation', 6 : 'infomap', 7 : 'infomap', 8 : 'modularity_optimization'}
##graph_file_separator = "\t"
scratch_directory_stem = "" #to be constructed after input parameters are processed
out_file_prefix = "clustering_v"
out_file_suffix = ".dat"
out_file_separator = '\t'
written_timing_data_header = False
timing_data_file_name = '' #to be constructed after input parameters are processed
logfile_name = '' #to be constructed after input parameters are processsed

####################
# Helper Functions #
####################

def handleArgs():
    """Handle command-line input arguments."""

    parser = argparse.ArgumentParser(description="Cluster a set of network graphs with various clustering algorithms.")
    parser.add_argument("-m", "--methods", nargs="+", type=str.lower, choices=supported_methods, required=True, help="the names of the clustering methods to run", dest="clustering_methods")
    parser.add_argument("-t", "--times", default=1, type=int, help="the number of clusterings to run on each graph", dest="numtrials")
    directionality_group = parser.add_mutually_exclusive_group(required=True)
    directionality_group.add_argument("-d", "--directed", action="store_true", help="indicates that the graphs are directed", dest="is_directed")
    directionality_group.add_argument("-u", "--undirected", action="store_true", help="indicates that the graphs are undirected", dest="is_undirected")
    parser.add_argument("--gpre", required=True, help="the stem for the path and filename of the graph files, before the file number", dest="graph_file_prefix")
    parser.add_argument("--gsuf", required=True, help="the ending to the filename of the graph files, after the file number, including the file extension", dest="graph_file_suffix")
    parser.add_argument("-s", "--start", default=1, type=int, help="the file number with which to start, defaults to 1", dest="file_range_start")
    parser.add_argument("-e", "--end", type=int, required=True, help="the file number with which to end, inclusive", dest="file_range_end")
    parser.add_argument("--lp", default=os.getcwd() + '/clustering_programs_5_2/', help="the path at which the Lancichinetti clustering program is installed, required if running the 'oslom', 'infomap', 'hierarchical_infomap', 'louvain', 'label_propagation', or 'modularity_optimization' methods. Defaults to the current working directory + '/clustering_programs_5_2/'", dest="lancichinetti_program_path")
    parser.add_argument("--modopt", default=os.getcwd(), help="the path at which Ludo Waltman and Nees Jan van Eck's Modularity Optimizer is installed, required if running the 'slm' method. Defaults to the current working directory", dest="modularity_optimizer_path")
    parser.add_argument("--Xmx", default="64m", help="A value for the Xmx Java parameter, to be passed to the call to the slm clustering algorithm. For example, if '600m' is given, the flag '-Xmx600m' will be passed to the java call to run the slm clustering algorithm. Defaults to '64m'", dest="xmx")
    parser.add_argument("--demon", default=os.getcwd() + '/demon_py', help="the path at which DEMON clustering, implemented by Giulio Rossetti, is installed, required if running the 'demon' method. Defaults to a folder called 'demon_py' in the current directory", dest="demon_path")
    parser.add_argument("-o", "--out", default="clustering_results/", help="the path to which to write the program output files, defaults to 'clustering_results/'", dest="out_directory_stem")
    
    global args
    args = parser.parse_args()

    args.xmx = '-Xmx' + args.xmx

def createPathIfNeeded(path):
    """Credits to user 'Heikki Toivonen' on SO: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary"""
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise

def deletePathIfNeeded(path):
    try:
        shutil.rmtree(scratch_directory_stem)
    except OSError as error:
        if error.errno != errno.ENOENT:
            raise

def deleteFileIfNeeded(file):
    try:
        os.remove(file)
    except OSError as error:
        if error.errno != errno.ENOENT:
            raise

def timing_values_from_err(err):

    timing_values = []

    timing_line = err.splitlines()[-2]
    assert ('user' in timing_line) and ('system' in timing_line) and ('elapsed' in timing_line)
    pieces = timing_line.split()

    user_piece = pieces[0] #i.e. user_piece = '0.46user'
    assert user_piece[-4:] == 'user'
    timing_values.append(user_piece[:-4])
    
    system_piece = pieces[1] #i.e. system_piece = '0.01system'
    assert system_piece[-6:] == 'system'
    timing_values.append(system_piece[:-6])

    elapsed_piece = pieces[2] #i.e. elapsed_piece = '0:02.51elapsed'
    assert elapsed_piece[-7:] == 'elapsed'
    timing_values.append(elapsed_piece[:-7])

    return timing_values

def write_timing_lines(method, network_number, trial_number, timing_values, write_header = True):

    write_to = args.out_directory_stem + timing_data_file_name

    #Writing to file here
    if write_header:
        with open(write_to, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Method', 'Network', 'Trial', 'User', 'System', 'Elapsed'])
        global written_timing_data_header
        written_timing_data_header = True

    with open(write_to, 'a') as f:
        writer = csv.writer(f)
        writer.writerow([method, network_number, trial_number, timing_values[0], timing_values[1], timing_values[2]])

def log_out_err(method, network_number, trial_number, out, err):

    write_to = args.out_directory_stem + logfile_name

    with open(write_to, 'a') as f:
        now = datetime.datetime.now().strftime('datetime: %X.%f on %x')
        f.write(now + '\n' + method + '_v' + network_number + '_' + trial_number + '\n\n' + out + '\n' + err + '\n\n')

def parseLancichinettiResults(f_path, out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path, p):

    # clustering file in (f_path + 'results_1/tp')

    read_file = open(f_path + 'results_1/tp', 'r')
    write_file = open(out_path + out_file_prefix + str(out_file_number) + out_file_suffix, 'wb')

    cluster_number_string = '1'
    for line in read_file:
        if line[0] == '#':
            continue
        nodes = line.split()
        for node in nodes:
            write_file.write(node + out_file_separator + cluster_number_string + '\n')
        cluster_number_string = str(int(cluster_number_string) + 1)

    print('\nSuccessfully ran ' + lancichinetti_num_to_name[p] + ' clustering and wrote results to file ' + out_path + out_file_prefix + str(out_file_number) + out_file_suffix + '\n')

    read_file.close()
    write_file.close()

def parseLeidenResults(leiden_file, min_node_id, out_file_prefix, out_file_number, out_file_suffix, out_file_separaotr, out_path):

    read_file = open(leiden_file, 'r')
    write_file = open(out_path + out_file_prefix + str(out_file_number) + out_file_suffix, 'wb')

    current_id = min_node_id
    for line in read_file:
        cluster_assignment = line.split()[0]
        write_file.write(str(current_id) + out_file_separator + cluster_assignment + '\n')
        current_id += 1

    print('\nSuccessfully ran Leiden clustering and wrote results to file ' + out_path + out_file_prefix + str(out_file_number) + out_file_suffix + '\n')

    read_file.close()
    write_file.close()

def parseDemonResults(demon_file, out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path):

    read_file = open(demon_file, 'r')
    write_file = open(out_path + out_file_prefix + str(out_file_number) + out_file_suffix, 'wb')

    node_set = set()
    for line in read_file:
        parts = line.split()
        cluster_num = str(int(parts[0])+ 1)
        nodes = parts[1].split(',')
        # Take only first assignment of node to cluster in file
        for node in nodes:
            if not node in node_set:
                node_set.add(node)
                write_file.write(node + out_file_separator + cluster_num + '\n')

    print('\nSuccessfully ran Demon clustering and wrote results to file ' + out_path + out_file_prefix + str(out_file_number) + out_file_suffix + '\n')

    read_file.close()
    write_file.close()

def clusterByLancichinetti(n, p, f, c, program_path):
    """Describe the function."""
    
    deletePathIfNeeded(f)
    process = subprocess.Popen(['time', 'python', 'select.py', '-n', n, '-p', str(p), '-f', f, '-c', str(c)], cwd = program_path, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = process.communicate()

    times = timing_values_from_err(err)
    write_timing_lines(lancichinetti_num_to_name[p], str(i + args.file_range_start), str(t + 1), times, write_header = not written_timing_data_header)
    log_out_err(lancichinetti_num_to_name[p], str(i + args.file_range_start), str(t + 1), out, err)

def clusterByLeiden(input_file, output_file, random_seed, modularity_function = 1, resolution_parameter = 1.0, optimization_algorithm = 3, n_random_starts = 10, n_iterations = 10, print_output = 0):
    
    deleteFileIfNeeded(args.modularity_optimizer_path + '/' + output_file)
    process = subprocess.Popen(['time', 'java', args.xmx, '-jar', 'ModularityOptimizer.jar', input_file, output_file, str(modularity_function), str(resolution_parameter), str(optimization_algorithm), str(n_random_starts), str(n_iterations), str(random_seed), str(print_output)], cwd=args.modularity_optimizer_path, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = process.communicate()

    times = timing_values_from_err(err)
    write_timing_lines('slm', str(i + args.file_range_start), str(t + 1), times, write_header = not written_timing_data_header)
    log_out_err('slm', str(i + args.file_range_start), str(t + 1), out, err)

def clusterByDemon(input_file):

    process = subprocess.Popen(['python', 'launch.py', input_file], cwd=args.demon_path, stdout = subprocess.PIPE, stderr = subprocess.PIPE)
    out, err = process.communicate()
    
    times = timing_values_from_err(err)
    write_timing_lines('demon', str(i + args.file_range_start), str(t + 1), times, write_header = not written_timing_data_header)
    log_out_err('demon', str(i + args.file_range_start), str(t + 1), out, err)

def runLancichinettiClustering(n, p, f, c, program_path, out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path):

    clusterByLancichinetti(n, p, f, c, program_path)
    parseLancichinettiResults(f, out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path, p)
    deletePathIfNeeded(scratch_directory_stem)

def runLeidenClustering(input_file, output_file, random_seed, out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path):

    clusterByLeiden(input_file, output_file, random_seed)
    parseLeidenResults(output_file, 0, out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path)
    deleteFileIfNeeded(args.modularity_optimizer_path + '/' + output_file)

def runDemonClustering(input_file, out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path):

    clusterByDemon(input_file) # Creates file demon_py/communities that needs to be parsed and that assigns some nodes to multiple communities
    parseDemonResults(args.demon_path + '/communities', out_file_prefix, out_file_number, out_file_suffix, out_file_separator, out_path)

if __name__ == "__main__":

    ##############################
    # Input Parameter Processing #
    ##############################

    handleArgs()

    scratch_directory_stem = 'scratch_folder_s_' + str(args.file_range_start) + '_e_' + str(args.file_range_end) + '/'
    timing_data_file_name = 'timing_data_s_' + str(args.file_range_start) + '_e_' + str(args.file_range_end) + '.csv'
    logfile_name = 'cluster_logfile_s_' + str(args.file_range_start) + '_e_' + str(args.file_range_end) + '.log'

    if args.is_directed:
        is_directed = True
    else:
        assert args.is_undirected
        is_directed = False

    graph_files = []
    for i in range(args.file_range_start, args.file_range_end + 1):
        graph_files.append(args.graph_file_prefix + str(i) + args.graph_file_suffix)

    ##################
    # Main Execution #
    ##################

    createPathIfNeeded(args.out_directory_stem)
    
    for method in args.clustering_methods:
        
        for i in xrange(len(graph_files)):

            for t in xrange(args.numtrials):

                if method == "blondel" or method == "louvain":
                    runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 4, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                elif method == "label_propagation":
                    runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 5, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                elif method == "modularity_optimization":
                    runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 8, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                elif method == "slm":
                    runLeidenClustering(os.getcwd() + '/' + graph_files[i], args.lancichinetti_program_path + 'modopt_output.txt', random.randrange(0, 1000000000000000000), method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)
                    
                elif is_directed:
                    if method == "oslom":
                        runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 1, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                    elif method == "infomap":
                        runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 3, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                    elif method == "hierarchical_infomap":
                        raise
                        #Implement correctly for hierarchical
                        runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 7, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)
                        
                else:
                    if method == "oslom":
                        runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 0, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                    elif method == "infomap":
                        runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 2, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                    elif method == "hierarchical_infomap":
                        raise
                        #Implement correctly for hierarchical
                        runLancichinettiClustering(os.getcwd() + '/' + graph_files[i], 6, os.getcwd() + '/' + scratch_directory_stem, 1, args.lancichinetti_program_path, method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)

                    elif method == "demon":
                        runDemonClustering(os.getcwd() + '/' + graph_files[i], method + "_" + out_file_prefix, str(i + args.file_range_start) + '_' + str(t + 1), out_file_suffix, out_file_separator, args.out_directory_stem)
