## Author: Scott Emmons (scott@scottemmons.com)
## Purpose: To increment Lancichinett's LFR graph generation program by the number of timeseed incrementation tokens that exist.
## Date: July 10, 2015

import argparse
import os
import glob

token_prefix = 'timeseed_inc_token'

def handleArgs():
    """Handles the command-line input arguments, placing them in the global Namespace variable 'args'."""

    parser = argparse.ArgumentParser(description="Increments Lancichinetti's LFR graph generation program by the number of timeseed incrementation tokens that exist")
    parser.add_argument("--wd", default=os.getcwd(), help="the working directory from which to count the number of timeseed incrementation tokens that exist", dest="working_directory")
    parser.add_argument("-b", "--benchmark", default=os.getcwd()+"/binary_networks/", help="the path to the installed LFR generation software to increment. Defaults to the current working directory + '/binary_networks/'", dest="bench_directory_stem")

    global args
    args = parser.parse_args()

def incrementBenchSeed(bench_directory_stem, times = 1):
    """Increments the value of the random seed in Lancichinetti's bencmark generation program."""

    seed_file = bench_directory_stem + 'time_seed.dat'
    
    with open(seed_file, 'r') as f:
        value = f.readline().split()[0]

    with open(seed_file, 'wb') as f:
        f.write(str(int(value) + 1 * times) + '\n')

handleArgs()

token_files = glob.glob(args.working_directory + '/' + token_prefix + '*')
if len(token_files) == 0:
    print 'WARNING: Zero token files were found.'    
incrementBenchSeed(args.bench_directory_stem, times = len(token_files))
for token_file in token_files:
    os.remove(token_file)
