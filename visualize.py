## Author: Scott Emmons
## Purpose: To visualize results of the clustering analysis workflow.
## Date: February 17, 2014

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt

####################
# Global variables #
####################

out_directory = 'generated_visualizations/'
compiled_data = out_directory + 'compiled_data.csv'

def handleArgs():
    """Handle command-line input arguments."""

    # Collect user input of raw data files to analyze, their corresponding names, and the number of nodes that corresponds to each file
    parser = argparse.ArgumentParser(description="Visualize the summary statistics of the clustering analysis workflow")
    parser.add_argument("-f", "--files", nargs="+", required=True,  help="the raw data files to visualize", dest="raw_data_file_paths")
    parser.add_argument("--names", nargs="+", required=True, help="the name of each raw data file, in corresponding order", dest="raw_data_file_names")
    parser.add_argument("--nodes", nargs="+", type=int, required=True, help="the number of nodes in each of the indicated raw data files, in corresponding order", dest="raw_data_nodes")

    global args
    args = parser.parse_args()

    # Ensure that the given raw data files, their names, and corresponding number of nodes are equal in length
    if len(args.raw_data_file_paths) != len(args.raw_data_file_names) != len(args.raw_data_nodes):
        print 'the length of -f (--files) must match that of --names and --nodes'
        assert False

def makeDirIfNeeded(directory):
    """"""

    try:
        os.mkdir(directory)
    except OSError:
        pass

def compileData(file_paths, destination_path, headers = True):
    """Write data from multiple data files into one file."""

    write_header = headers

    makeDirIfNeeded(out_directory)
    destination = open(destination_path, 'wb')

    for path in file_paths:
        source = open(path, 'r')
        # Skip column names for all but first data file
        if not write_header:
            source.readline()
        else:
            write_header = False
        for line in source:
            destination.write(line)
        source.close()

    destination.close()

    status_message = '\nSuccessfully compiled data from file(s)'
    for path in file_paths:
        status_message = status_message + ' ' + path + ','
    print status_message[:-1]
    
def runIndividualVisualization(file_path, file_name):
    """
    Visualize the individual metrics in file_path, using file_name to name the output.
    Dependencies:
    pandas Python library
    Input variables:
    file_path - a string specifying the path to the file containing the data to visualize
    file_name - a string specifying the name of the file to use for function output naming
    Output:
    writes a file visualizing the metrics contained in file_path
    """

    # Create a data frame of the data stored in the file
    df = pd.read_csv(file_path)

    # For each clustering algorithm, visualize performance by metric
    algorithms = df['name'].unique()
    for algorithm in algorithms:
        algorithm_frame = df[df['name'] == algorithm][['metric', 'value']]
        boxplot = algorithm_frame.boxplot(by="metric", showmeans=True, meanline=True)
        boxplot.set_title(algorithm + ' Clustering Metric Values')
        boxplot.set_xlabel('')
        boxplot.tick_params(axis="x", labelsize=5)
        boxplot.set_ylim(-0.1, 1.1)
        figure = boxplot.get_figure()
        figure.suptitle(file_name)
        makeDirIfNeeded(out_directory + file_name + '/')
        figure.savefig(out_directory + file_name + '/' + algorithm + '_boxplot' + '.pdf')
        plt.close()

    # For each evaluation metric, visualize performance by clustering algorithm
    metrics = df['metric'].unique()
    for metric in metrics:
        metric_frame = df[df['metric'] == metric][['name', 'value']]
        boxplot = metric_frame.boxplot(by="name", showmeans=True, meanline=True)
        boxplot.set_title(metric + ' Performance by Algorithm')
        boxplot.set_xlabel('')
        boxplot.set_ylim(-0.1, 1.1)
        figure = boxplot.get_figure()
        figure.suptitle(file_name)
        makeDirIfNeeded(out_directory + file_name + '/')
        figure.savefig(out_directory + file_name + '/' + metric.lower().replace(' ', '_') + '_boxplot' + '.pdf')
        plt.close()

    print('\nSuccessfully ran individual visualization on data file ' + file_path)

def runCombinedVisualization(combined_data, write_subdirectory = ''):
    """"""

    # Create a data frame of the data stored in the file
    df = pd.read_csv(combined_data)

    # For each clustering algorithm, visualize performance by metric
    algorithms = df['name'].unique()
    for algorithm in algorithms:
        algorithm_frame = df[df['name'] == algorithm][['n', 'metric', 'value']]
        
        # Create dictionary to pass to dataframe constructor for visualization
        d = {}

        for metric in algorithm_frame['metric'].unique():
            series = pd.Series()
            metric_match = algorithm_frame['metric'] == metric
            for n in algorithm_frame['n'].unique():
                n_match = algorithm_frame['n'] == n
                dat = algorithm_frame[metric_match & n_match]['value']
                series = series.append(pd.Series([dat.mean()], index=[n]))
            d[metric] = series
        
        # Create and visualize the dataframe
        vis_frame = pd.DataFrame(d)

        linegraph = vis_frame.plot(marker='^')
        linegraph.set_title(algorithm + ' Clustering Metric Values')
        linegraph.set_xscale('log')
        linegraph.set_xlim(linegraph.get_xlim()[0] * 0.90, linegraph.get_xlim()[1] * 1.1)
        linegraph.set_ylim(linegraph.get_ylim()[1] * -0.10, linegraph.get_ylim()[1] * 1.1)
        figure = linegraph.get_figure()
        figure.tight_layout()
        makeDirIfNeeded(out_directory + write_subdirectory)
        figure.savefig(out_directory + write_subdirectory + '/' + algorithm.lower().replace(' ', '_') + '_linegraph' + '.pdf')
        plt.close()

    # For each evaluation metric, visualize performance by clustering algorithm
    metrics = df['metric'].unique()
    for metric in metrics:
        metric_frame = df[df['metric'] == metric][['n', 'name', 'value']]
        
        # Create dictionary to pass to dataframe constructor for visualization
        d = {}

        for algorithm in metric_frame['name'].unique():
            series = pd.Series()
            algorithm_match = metric_frame['name'] == algorithm
            for n in metric_frame['n'].unique():
                n_match = metric_frame['n'] == n
                dat = metric_frame[algorithm_match & n_match]['value']
                series = series.append(pd.Series([dat.mean()], index=[n]))
            d[algorithm] = series
        
        # Create and visualize the dataframe
        vis_frame = pd.DataFrame(d)

        linegraph = vis_frame.plot(marker='^')
        linegraph.set_title(metric + ' Performance by Algorithm')
        linegraph.set_xscale('log')
        linegraph.set_xlim(linegraph.get_xlim()[0] * 0.90, linegraph.get_xlim()[1] * 1.1)
        linegraph.set_ylim(-0.1, 1.1)
        figure = linegraph.get_figure()
        figure.tight_layout()
        makeDirIfNeeded(out_directory + write_subdirectory)
        figure.savefig(out_directory + write_subdirectory + '/' + metric.lower().replace(' ', '_') + '_linegraph' + '.pdf')
        plt.close()

    print('\nSuccessfully ran combined visualization on data file ' + combined_data)

if __name__ == "__main__":

    # Parse command-line arguments
    handleArgs()
    
    # Compile raw data file paths into single file for integrated analysis
    compileData(args.raw_data_file_paths, compiled_data)

    # Iterate over each raw data file and create individual visualizations
    for i in xrange(len(args.raw_data_file_paths)):

        file_path = args.raw_data_file_paths[i]
        file_name = args.raw_data_file_names[i]
        nodes = args.raw_data_nodes[i]
        runIndividualVisualization(file_path, file_name)

    # Visualize combined data to show how algorithms scale with n
    runCombinedVisualization(compiled_data, write_subdirectory = 'combined')
