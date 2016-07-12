## Author: Scott Emmons
## Purpose: To visualize results of the clustering analysis workflow.
## Date: February 17, 2015

import argparse
import os
import pandas as pd
import numpy as np
import matplotlib
matplotlib.use('Agg')
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
    

def formatName(name):
    """"""

    if name == 'blondel':
	return 'Louvain'
    elif name == 'infomap':
	return 'Infomap'
    elif name == 'label_propagation':
	return 'Label Propagation'
    elif name == 'slm':
	return 'Smart Local Moving'
    elif name == 'Scikit-learn NMI':
	return 'Traditional NMI'
    elif name == 'Lancichinetti NMI':
	return 'Variant of NMI'
    
    return name

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
	mean_map = {}
	error_map = {}

        for metric in algorithm_frame['metric'].unique():
            means = pd.Series()
            errors = pd.Series()
            metric_match = algorithm_frame['metric'] == metric
            means = means.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
            errors = errors.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
            for n in algorithm_frame['n'].unique():
                n_match = algorithm_frame['n'] == n
                dat = algorithm_frame[metric_match & n_match]['value']
                means = means.append(pd.Series([dat.mean()], index=[n]))
                errors = errors.append(pd.Series([dat.std()], index=[n]))
            mean_map[metric] = means
            error_map[metric] = errors
        
        # Create dataframes of relevant information
        mean_frame = pd.DataFrame(mean_map)
        error_frame = pd.DataFrame(error_map)

	# Create line graph with standard deviations
        linegraph = mean_frame.plot(marker='^', yerr=error_frame)
        linegraph.set_title(algorithm + ' Clustering Metric Values')
        linegraph.set_xscale('log')
        linegraph.set_xlim(sorted(mean_frame.index.values)[1] * 0.90, linegraph.get_xlim()[1] * 1.1) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
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
        mean_map = {}
	error_map = {}

        for algorithm in metric_frame['name'].unique():
            means = pd.Series()
            errors = pd.Series()
            algorithm_match = metric_frame['name'] == algorithm
            means = means.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
            errors = errors.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
            for n in metric_frame['n'].unique():
                n_match = metric_frame['n'] == n
                dat = metric_frame[algorithm_match & n_match]['value']
                means = means.append(pd.Series([dat.mean()], index=[n]))
                errors = errors.append(pd.Series([dat.std()], index=[n]))
            mean_map[algorithm] = means
            error_map[algorithm] = errors
        
        # Create dataframes of relevant information
        mean_frame = pd.DataFrame(mean_map)
        error_frame = pd.DataFrame(error_map)

	# Create line graph with standard deviations
        linegraph = mean_frame.plot(marker='^', yerr=error_frame)
        linegraph.set_title(metric + ' Performance by Algorithm')
        linegraph.set_xscale('log')
        linegraph.set_xlim(sorted(mean_frame.index.values)[1] * 0.90, linegraph.get_xlim()[1] * 1.1) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
        linegraph.set_ylim(-0.1, 1.1)
        figure = linegraph.get_figure()
        figure.tight_layout()
        makeDirIfNeeded(out_directory + write_subdirectory)
        figure.savefig(out_directory + write_subdirectory + '/' + metric.lower().replace(' ', '_') + '_linegraph' + '.pdf')
        plt.close()

    # Create algorithm by metric subplot matrix
    algorithms = [algorithm for algorithm in algorithms if algorithm != 'gold_standard'] 
    fig, axes = plt.subplots(len(algorithms), len(metrics), sharex='col', sharey='row')
    fig.suptitle(r"LFR Benchmark Results")# at $\mu$ = 0.40")
    for i in range(len(algorithms)):
	for j in range(len(metrics)):
	    paired_frame = df[(df['name'] == algorithms[i]) & (df['metric'] == metrics[j])][['n', 'value']]
	    dataset = [paired_frame[paired_frame['n'] == n][['value']].values for n in df['n'].unique()]
	    axes[i, j].violinplot(dataset, showmeans=True, showextrema=True)
	    axes[i, j].set_ylim(-0.1, 1.1)
	    for tick in axes[i, j].get_xticklabels():
	        tick.set_rotation(90)
		tick.set_fontsize(8)
	    if (i == 0):
		axes[i, j].set_title(formatName(metrics[j]), fontsize=7)
	    if (j == 0):
		axes[i, j].set_ylabel(formatName(algorithms[i]), rotation=90, fontsize=9)
    plt.setp(axes, xticks=[x + 1 for x in range(df['n'].unique().size)], xticklabels=df['n'].unique().tolist())
    makeDirIfNeeded(out_directory + write_subdirectory)
    fig.savefig(out_directory + write_subdirectory + '/' + 'violinplot_matrix' + '.pdf')
    plt.close()

    print('\nSuccessfully ran combined visualization on data file ' + combined_data)

def runSummaryVisualization(combined_data, write_subdirectory=''):
    """"""

    # Create a data frame of the data stored in the file
    df = pd.read_csv(combined_data)
    df = df[df['name'] != 'gold_standard']

    # Visualize performance of algorithms averaged over choice of metric
    algorithm_frame = df[['n', 'name', 'value']]

    # Create dictionary to pass to dataframe constructor for visualization
    mean_map = {}
    error_map = {}
    
    for algorithm in algorithm_frame['name'].unique():
        means = pd.Series()
        errors = pd.Series()
        algorithm_match = algorithm_frame['name'] == algorithm
        means = means.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
        errors = errors.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
        for n in algorithm_frame['n'].unique():
            n_match = algorithm_frame['n'] == n
            dat = algorithm_frame[algorithm_match & n_match]['value']
            means = means.append(pd.Series([dat.mean()], index=[n]))
            errors = errors.append(pd.Series([dat.std()], index=[n]))
        mean_map[algorithm] = means
        error_map[algorithm] = errors

    # Create dataframes of relevant information
    mean_frame = pd.DataFrame(mean_map)
    error_frame = pd.DataFrame(error_map)

    # Create line graph with standard deviations
    linegraph = mean_frame.plot(marker='^', yerr=error_frame)
    linegraph.set_title('Algorithm Performance Averaged Over All Metrics')
    linegraph.set_xscale('log')
    linegraph.set_xlim(sorted(mean_frame.index.values)[1] * 0.90, linegraph.get_xlim()[1] * 1.1) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
    linegraph.set_ylim(-0.1, 1.1)
    figure = linegraph.get_figure()
    figure.tight_layout()
    makeDirIfNeeded(out_directory + write_subdirectory)
    figure.savefig(out_directory + write_subdirectory + '/' + 'algorithm_linegraph' + '.pdf')
    plt.close()

    # Visualize performance by metrics average over choice of algorithm
    metric_frame = df[['n', 'metric', 'value']]

    # Create dictionary to pass to dataframe constructor for visualization
    mean_map = {}
    error_map = {}

    for metric in metric_frame['metric'].unique():
        means = pd.Series()
        errors = pd.Series()
        metric_match = metric_frame['metric'] == metric
        means = means.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
        errors = errors.append(pd.Series([0], index=[0])) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
        for n in metric_frame['n'].unique():
            n_match = metric_frame['n'] == n
            dat = metric_frame[metric_match & n_match]['value']
            means = means.append(pd.Series([dat.mean()], index=[n]))
            errors = errors.append(pd.Series([dat.std()], index=[n]))
        mean_map[metric] = means
        error_map[metric] = errors

    # Create dataframes of relevant information
    mean_frame = pd.DataFrame(mean_map)
    error_frame = pd.DataFrame(error_map)

    # Create line graph with standard deviations
    linegraph = mean_frame.plot(marker='^', yerr=error_frame)
    linegraph.set_title('Metric Performance Averaged Over All Algorithms')
    linegraph.set_xscale('log')
    linegraph.set_xlim(sorted(mean_frame.index.values)[1] * 0.90, linegraph.get_xlim()[1] * 1.1) #Dirty solution to side step bug https://github.com/pydata/pandas/issues/11858
    linegraph.set_ylim(-0.1, 1.1)
    figure = linegraph.get_figure()
    figure.tight_layout()
    makeDirIfNeeded(out_directory + write_subdirectory)
    figure.savefig(out_directory + write_subdirectory + '/' + 'metric_linegraph' + '.pdf')
    plt.close()

    print('\nSuccessfully ran summary visualization on data file ' + combined_data)

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

    # Visualize summary graphs to picture overall trends
    runSummaryVisualization(compiled_data, write_subdirectory = 'summary')
