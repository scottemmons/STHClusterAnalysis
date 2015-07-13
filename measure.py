## Author: Scott Emmons
## Purpose: To evaluate clustering metrics on the input network graph and partitions.
## Date: January 7, 2015

## Instructions:
## This program will evaluate clustering metrics over multiple pairs of graphs and partitions.
## Each input graph file should be an edgelist of the form "source" + "separator" + "destination"
## Each input partition file should be of the form "node" + "separator" + "assignment"
## All of the input graph files and input partition should be named with a consistent convention such that the naming scheme for each type of file is consistent.
## The only variation in the names of the input files should be a number within an integral range that pairs the files.
## This program will evaluate metrics for each pair of graph and partition files within the given range.

## Parameters:
## graph_file_prefix is a string specifying the naming convention for graph files before the file identificatio number
## graph_file_suffix is a string specifying the naming convention for graph files after the file identification number
## graph_file_separator is a string serving as the graph file column separator
## partition_file_prefix is a string specifying the naming convention for partition files before the file identificatio number
## partition_file_suffix is a string specifying the naming convention for partition files after the file identification number
## partition_file_separator is a string serving as the partition file column separator
## file_range_start marks the beginning of the integral identification number range, inclusive
## file_range_end marks the end of the integral identification number range, inclusive
## output_path is a string specifying the path to which the output files will be written
## output_file_prefix is a string specifying the naming convention for the output files before the file identification number
## output_file_suffix is a string specifying the naming convention for the output files before the file identification number

## Notes to self:
## Lancichinetti NMI for covers: https://sites.google.com/site/andrealancichinetti/mutual
## GMap metric definitions paper: https://www.cs.arizona.edu/~kobourov/contiguous.pdf

## Testing command: python measure.py --gmap gmap/ --gpre real_world_graphs/flickr/flickrEdges_renum_v --gsuf .txt --gsep $'\t' -u -s 1 -e 1 --cnames blondel infomap label_propagation slm --cpre real_world_graphs/flickr/blondel_clustering_v real_world_graphs/flickr/infomap_clustering_v real_world_graphs/flickr/label_propagation_clustering_v real_world_graphs/flickr/slm_clustering_v --csuf .dat --csep $'\t' --cnum 1 -o real_world_graphs/flickr/

import argparse
import os
import errno
import math
import csv
import subprocess

import igraph
from sklearn.metrics.cluster import normalized_mutual_info_score
from sklearn.metrics.cluster import adjusted_rand_score

####################
# Global Variables #
####################

file_range_start_default = 1
output_path_default = "metric_results/"
output_file_prefix = 'results_v'
output_file_suffix = '.csv'
dot_file_prefix = 'dotfile_v'
dot_file_suffix = '.dot'
raw_data_file_name = '' #to be constructed after input parameters are processed
written_raw_data_header = False
summary_file_name = 'summary_statistics.csv'
logfile_lines = []

####################
# Helper functions #
####################

def handleArgs():
    """Handle command-line input arguments."""

    parser = argparse.ArgumentParser(description="Run metrics of the quality of graph clusterings over a set of network graphs and clustering partitions.")
    parser.add_argument("--gmap", required=True, help="path to an installed version of GMap's external metric analysis program, such as can be found at https://github.com/spupyrev/gmap/tree/master/external/eba", dest="gmap_directory")
    parser.add_argument("--lnmi", default="", help="path to an installed version of Andrea Lancichinetti's normalized mutual information software, which can be found at https://sites.google.com/site/andrealancichinetti/mutual", dest="lnmi_directory")
    parser.add_argument("--gpre", required=True, help="the stem for the path and filename of the graph files, before the file number", dest="graph_file_prefix")
    parser.add_argument("--gsuf", required=True, help="the ending to the filename of the graph files, after the file number, including the file extension", dest="graph_file_suffix")
    parser.add_argument("--gsep", required=True, help="the column separator in the graph files", dest="graph_file_separator")
    directionality_group = parser.add_mutually_exclusive_group(required=True)
    directionality_group.add_argument("-d", "--directed", action="store_true", help="indicates that the graphs are directed", dest="is_directed")
    directionality_group.add_argument("-u", "--undirected", action="store_true", help="indicates that the graphs are undirected", dest="is_undirected")
    parser.add_argument("--srun", action="store_true", help="indicates to run gold standard analysis", dest="run_standard")
    parser.add_argument("--spre", default="", help="the stem for the path and filename of the 'gold_standard' clutsering files, before the file number", dest="gold_standard_file_prefix")
    parser.add_argument("--ssuf", default="", help="the ending to the filename of the 'gold standard' clustering files, after the file number, including the file extension", dest="gold_standard_file_suffix")
    parser.add_argument("--ssep", default="", help="the column separator in the 'gold standard' clutsering files", dest="gold_standard_file_separator")
    parser.add_argument("-s", "--start", default=file_range_start_default, type=int, help="the file number with which to start, defaults to 1", dest="file_range_start")
    parser.add_argument("-e", "--end", type=int, required=True, help="the file number with which to end, inclusive", dest="file_range_end")
    parser.add_argument("--cnames", nargs="+", default=[], type=str.lower, required=True, help="the names of the clustering methods that will be evaluated, to be used in the naming of output files", dest="clustering_file_names")
    parser.add_argument("--cpre", nargs="+", default=[], help="the stem for the path and filename of the to-be-evaluated clustering files, before the file number", dest="clustering_file_prefixes")
    parser.add_argument("--csuf", nargs="+", default=[], help="the ending to the filename of the to-be-evaluated clustering files, after the file number, including the file extension; must be either a list matching the length of --cpre, or one value that is universal to all in --cpre", dest="clustering_file_suffixes")
    parser.add_argument("--csep", nargs="+", default=[], help="the column separator in the to-be-evaluated clustering files; must be either a list matching the length of --cpre, or one value that is universal to all in --cpre", dest="clustering_file_separators")
    parser.add_argument("--cnum", type=int, default=1, help="the number of clusterings that exist for each graph", dest="clusterings_per_graph")
    parser.add_argument("-o", "--out", default=output_path_default, help="the directory to which to write the program output files, defaults to 'metric_results/'", dest="output_path")

    global args
    args = parser.parse_args()

    if len(args.clustering_file_names) != len(args.clustering_file_prefixes):
        print 'the length of --cnames must match that of --cpre'
        assert False

    if len(args.clustering_file_prefixes) != len(args.clustering_file_suffixes):
        if len(args.clustering_file_suffixes) != 1:
            print 'the length of --csuf must either match that of --cpre or be equal to one, specifying a universal suffix'
            assert False

    if len(args.clustering_file_prefixes) != len(args.clustering_file_separators):
        if len(args.clustering_file_separators) != 1:
            print 'the length of --csep must either match that of --cpre or be equal to one, specifying a universal separator'
            assert False

    if args.run_standard:
        try:
            assert (len(args.lnmi_directory) > 0) and (len(args.gold_standard_file_prefix) > 0) and (len(args.gold_standard_file_suffix) > 0) and (len(args.gold_standard_file_separator) > 0)
        except:
            print "If --srun flag is given, the program will attempt to run gold standard metric analysis and you must supply --lnmi, --spre, --ssuf, and --ssep"
            assert False

def writeLinesToFile(lines, filename, mode = 'wb'):
    """"""

    with open(filename, mode) as f:
        for line in lines:
            f.write(line + '\n')

def appendLines(to_add, add_to):
    """Append lines to_add, a list, to add_to, a list."""
    
    for line in to_add:
        add_to.append(line)

def createPathIfNeeded(path):
    """Credits to user 'Heikki Toivonen' on SO: http://stackoverflow.com/questions/273192/check-if-a-directory-exists-and-create-it-if-necessary"""
    try:
        os.makedirs(path)
    except OSError as error:
        if error.errno != errno.EEXIST:
            raise

def partitionFromFile(partition_file, partition_file_separator):
    """Create a partition object from the given file.
    Return dictionary object assigning nodes to clusters based on
    partition_file, a string.
    Returns dictionary object assigning nodes to clusters."""

    partition = {}

    f = open(partition_file, 'r')
    for line in f:
        node, cluster = line.split(partition_file_separator)
        partition[int(node)] = int(cluster.rstrip())
    f.close()

    return partition

def writeDotFile(graph_file_name, graph_separator, partition_file_name, partition_separator, output_path, is_directed = True):
    """"""

    write_file = open(output_path, 'wb')
    write_file.write('graph {\n')

    # Write node id lines
    with open(partition_file_name, 'r') as partition_file:
        for line in partition_file:
            pieces = line.split(partition_separator)
            id = pieces[0]
            cluster = pieces[1].rstrip() # remove newline character and trailing spaces
            write_file.write('  "' + id + '" [cluster="' + cluster + '"];\n')

    # Write edge lines
    graph_file = open(graph_file_name, 'r')
##    if is_directed:
##        edge_str = '->'
##    else:
##        edge_str = '--'
    edge_str = '--'
    redundant_edges = {}
    for line in graph_file:
        pieces = line.split(graph_separator)
        source = pieces[0]
        destination = pieces[1].rstrip() # remove newline character and trailing spaces
        try:
            redundant_edges[destination].add(source)
        except KeyError:
            redundant_edges[destination] = set()
            redundant_edges[destination].add(source)
        if not destination in redundant_edges.get(source, set()):
            write_file.write('  "' + source + '" ' + edge_str + ' "' + destination + '";\n')
    graph_file.close()

    write_file.write('}')
    write_file.close()

def mergeOccurrenceDicts(dict_1, dict_2):
    """"""

    merged_dict = {}

    unique_keys = set()
    for key in dict_1.keys():
        unique_keys.add(key)
    for key in dict_2.keys():
        unique_keys.add(key)

    for key in unique_keys:
        merged_dict[key] = dict_1.get(key, 0) + dict_2.get(key, 0)

    return merged_dict

def writeRawDataLines(lines, name, n, network_id, clustering_id):
    """"""

    global written_raw_data_header
    write_to = args.output_path + raw_data_file_name

    if not written_raw_data_header:
        with open(write_to, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['name', 'metric', 'cluster_num', 'value', 'n', 'network_id', 'clustering_id'])
            written_raw_data_header = True

    with open(write_to, 'a') as f:
        writer = csv.writer(f)
        for line in lines:
            writer.writerow([name] + line + [n, network_id, clustering_id])

    print '\nResults for ' + name + ' on iteration ' + str(network_id) + '_' + str(clustering_id) + ' appended to ' + os.getcwd() + '/' + write_to

def exportSummaryLines(lines, name, write_header = True):
    """"""

    write_to = args.output_path + summary_file_name

    if write_header:
        with open(write_to, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Name', 'Metric', 'Statistic', 'Value'])

    with open(write_to, 'a') as f:
        writer = csv.writer(f)
        for line in lines:
            writer.writerow(line)

    print '\nSuccessfully wrote summary lines for ' + name + ' to ' + os.getcwd() + '/' + write_to

def calcOccurrenceDictMedian(dict, datasize = None):
    """"""

    values = dict.keys()
    values.sort()

    if not datasize:
        datasize = 0.0
        for value in values:
            datasize += value * dict[value]
    else:
        float(datasize)

    count = 0.0
    for i in xrange(len(values)):
        value = values[i]
        count += dict[value]
        difference = (datasize / 2.0) - count
        if difference < 0:
            median = value
            break
        elif difference == 0:
            median = (value + values[i+1]) / 2.0
            break

    return median

def calcOccurrenceDictStandardDeviation(dict, mean):
    """"""

    values = dict.keys()
    datasize = 0.0
    sigma_term = 0.0
    for value in values:
        sigma_term += float(dict[value]) * (value - mean)**2
        datasize += float(dict[value])
    return (sigma_term / datasize) ** (0.5) 

def statisticsFromOccurrenceDict(dict):
    """"""

    lines = []

    datasize = 0.0
    values = dict.keys()

    min = max = mean_term = float(values[0])
    datasize += dict[values[0]]
    mean_term *= datasize
    for i in xrange(1, len(values)):
        value = values[i]
        if value < min:
            min = value
        elif value > max:
            max = value
        mean_term += value * dict[value]
        datasize += dict[value]

    mean = mean_term / datasize
    median = calcOccurrenceDictMedian(dict, datasize = datasize)
    sd = calcOccurrenceDictStandardDeviation(dict, mean)

    lines.append(['minimum', min])
    lines.append(['maximum', max])
    lines.append(['median', median])
    lines.append(['mean', mean])
    lines.append(['standard deviation', sd])

    return lines

def calcVectorStandardDeviation(vector, mean):
    sigma_term = 0.0
    for value in vector:
        sigma_term += (value - mean)**2
    return (sigma_term / float(len(vector)))**(0.5)

def statisticsFromVector(values):
    """"""

    lines = []

    length = len(values)
    min = max = mean_term = float(values[0])
    for i in xrange(1, length): 
        value = values[i]
        mean_term += value
        if value < min:
            min = value
        elif value > max:
            max = value
    mean = mean_term / float(length)
    median = (values[int(math.floor(float(length - 1) / 2.0))] + values[int(math.ceil(float(length - 1) / 2.0))]) / 2.0
    sd = calcVectorStandardDeviation(values, mean)
    
    lines.append(['minimum', min])
    lines.append(['maximum', max])
    lines.append(['median', median])
    lines.append(['mean', mean])
    lines.append(['standard deviation', sd])

    return lines

def getStatLines(measurement_lookup, name):
    """"""

    lines = []

    for metric in measurement_lookup.keys():
        for line in statisticsFromVector(measurement_lookup[metric]):
            statistic, value = line
            lines.append([name, metric, statistic, value])

    return lines

def exportSummaryStatistics(data_lines, write_header = False):
    """"""

    write_to = args.output_path + summary_file_name
    if write_header:
        with open(write_to, 'w') as f:
            writer = csv.writer(f)
            writer.writerow(['Name', 'Metric', 'Statistic', 'Value'])

    for data in data_lines:

        name = data[0]
        measurements = {}
        for i in xrange(1, len(data)):
            trial = data[i]
            for measurement in trial:
                metric = measurement[0].lower()
                try:
                    value = float(measurement[2])
                except ValueError:
                    print 'Undefined value for metric:', metric
                    value = 0.0
                try:
                    measurements[metric].append(value)
                except KeyError:
                    measurements[metric] = [value]
        
        stat_lines = getStatLines(measurements, name)
        with open(write_to, 'a') as f:
            writer = csv.writer(f)
            writer.writerows(stat_lines)
        print '\nSuccessfully ran summary statistics for ' + name + ' and appended results to ' + os.getcwd() + '/' + write_to
        
def generateCorrespondingVectors(partition_1, partition_2):
    """"""

    vector_1 = []
    vector_2 = []

    for key in partition_1.keys():
        vector_1.append(partition_1[key])
        vector_2.append(partition_2[key])

    return vector_1, vector_2

def writeLineDefinedClusterFile(vector, to_write_path):
    """From vector of cluster definitions which is of the form
    node i assigned to value of vector[i], write the clutser assignments
    to a file assigning the nodes on a line to the same cluster."""

    cluster_to_node = {}

    for i in xrange(len(vector)):
        try:
            cluster_to_node[vector[i]].append(i + 1)
        except KeyError:
            cluster_to_node[vector[i]] = [i + 1]

    with open(to_write_path, 'wb') as f:
        for cluster in cluster_to_node.keys():
            node_list = cluster_to_node[cluster]
            first_node = True
            for node in node_list:
                if first_node:
                    f.write(str(node))
                    first_node = False
                else:
                    f.write(' ' + str(node))
            f.write('\n')

def runScikitNormMutInf(gold_standard_vector, partition_vector):
    """Normalized mutual information as defined here: http://scikit-learn.org/stable/modules/generated/sklearn.metrics.normalized_mutual_info_score.html#sklearn.metrics.normalized_mutual_info_score""" 

    result_lines = []

    value = normalized_mutual_info_score(gold_standard_vector, partition_vector)

    result_lines.append(['Scikit-learn NMI', 'Entire Graph', value])

    return result_lines

def runAdjRandScr(gold_standard_vector, partition_vector):
    """Adjusted Rand Score as defined here: http://scikit-learn.org/stable/modules/generated/sklearn.metrics.adjusted_rand_score.html#sklearn.metrics.adjusted_rand_score"""

    result_lines = []

    value = adjusted_rand_score(gold_standard_vector, partition_vector)

    result_lines.append(['Adjusted Rand Score', 'Entire Graph', value])

    return result_lines

def runLancichNormMutInf(gold_standard_vector, partition_vector):
    """Normalized Mutual Information as defined here: https://sites.google.com/site/andrealancichinetti/mutual"""

    result_lines = []

    gold_standard_cluster_name = "file1.dat"
    partition_cluster_name = "file2.dat"

    writeLineDefinedClusterFile(gold_standard_vector, args.lnmi_directory + gold_standard_cluster_name)
    writeLineDefinedClusterFile(partition_vector, args.lnmi_directory + partition_cluster_name)

    process = subprocess.Popen(['./mutual', gold_standard_cluster_name, partition_cluster_name], cwd=args.lnmi_directory, stdout=subprocess.PIPE)
    output = process.communicate()[0].split()
    assert output[0] == 'mutual3:'
    value = float(output[1])

    os.remove(args.lnmi_directory + gold_standard_cluster_name)
    os.remove(args.lnmi_directory + partition_cluster_name)

    result_lines.append(['Lancichinetti NMI', 'Entire Graph', value])

    return result_lines

def runGMapAnalysis(relative_dotfile_path):
    """"""

    result_lines = []

    modularity = None
    conductance = None
    coverage = None

    gmap_metric_directory = args.gmap_directory + 'external/eba'
    absolute_dotfile_path = os.getcwd() + '/' + relative_dotfile_path
    subprocess.call(['./kmeans', '-action=metrics', '-o=metric_results.txt', absolute_dotfile_path], cwd = gmap_metric_directory)

    gmap_output_path = gmap_metric_directory + '/metric_results.txt'
    f = open(gmap_output_path, 'r')
    for line in f:
        pieces = line.split()
        metric = pieces[0]
        value = pieces[1]
        if metric[:10] == 'Modularity':
            result_lines.append(['Modularity', 'Entire Graph', value])
        elif metric[:11] == 'Conductance':
            try:
                float(value)
                result_lines.append(['Conductance', 'Entire Graph', value])
            except ValueError:
                assert value == 'undefined'
                logfile_lines.append('Undefined conductance for ' + relative_dotfile_path)
                result_lines.append(['Conductance', 'Entire Graph', 0.0])
        elif metric[:8] == 'Coverage':
            result_lines.append(['Coverage', 'Entire Graph', value])

    f.close()
    os.remove(gmap_output_path)

    return result_lines

def runShortestPaths(graph):
    """Graph is an iGraph graph object."""

    occurrence_dict = {}

    shortest_paths = graph.shortest_paths()
    for i in range(len(shortest_paths)):
        matrix_row = shortest_paths[1]
        for j in range(len(matrix_row)):
            distance = matrix_row[j]
            if i != j and not math.isinf(distance) and distance != 0:
                assert type(distance) == int
                occurrence_dict[distance] = occurrence_dict.get(distance, 0) + 1

    return occurrence_dict

def runClusteringCoeff(graph):
    """Graph is an iGraph graph object."""

    return graph.transitivity_avglocal_undirected()

def runIndividualGraphMetrics(graph, partition, dotfile_path, is_directed):
    """"""

    result_lines = []
                        
##    result_lines = [["Metric", "Cluster Number", "Value"]]

    gmap_lines = runGMapAnalysis(dotfile_path)
    appendLines(gmap_lines, result_lines)

    return result_lines

def runComparisonGraphMetrics(graph, partition, gold_standard, is_directed):
    """"""

    result_lines = []


    gold_standard_vector, partition_vector = generateCorrespondingVectors(gold_standard, partition)

    # Metrics from scikit-learn

    scikit_norm_mut_inf_lines = runScikitNormMutInf(gold_standard_vector, partition_vector)
    appendLines(scikit_norm_mut_inf_lines, result_lines)

    adj_rand_scr_lines = runAdjRandScr(gold_standard_vector, partition_vector)
    appendLines(adj_rand_scr_lines, result_lines)

    # Lancichinetti's NMI measure

    lancich_norm_mut_inf_lines = runLancichNormMutInf(gold_standard_vector, partition_vector)
    appendLines(lancich_norm_mut_inf_lines, result_lines)

    return result_lines

def structuralMetricAnalysis(graph_file_list, graph_separator, name):
    """"""

    results = []

    shortest_paths_occurrence_dict = {}
    clustering_coeff_values = []

    for i in xrange(len(graph_file_list)):

        graph = igraph.Graph.Read_Edgelist(graph_file_list[i], directed = False)

        # Collect shortest path and clustering coefficient values
        shortest_paths_occurrence_dict = mergeOccurrenceDicts(runShortestPaths(graph), shortest_paths_occurrence_dict)
        clustering_coeff_values.append(runClusteringCoeff(graph))

        print '\nSuccessfully analyzed structural metrics of network for file number ' + str(args.file_range_start + i)

    # Parse collected values into result lines
    for line in statisticsFromOccurrenceDict(shortest_paths_occurrence_dict):
        statistic, value = line
        results.append([name, 'shortest_path', statistic, value])

    for line in statisticsFromVector(clustering_coeff_values):
        statistic, value = line
        results.append([name, 'clustering_coefficient', statistic, value])

    return results

def individualMetricAnalysis(graph_file_list, partition_file_lists, graph_separator, partition_separator, name):
    """"""
    
    results = [] 

    for i in xrange(((len(graph_file_list) + len(partition_file_lists)) / 2)): # assume graph_file_list and partition_file_lists are equal in length

        graph = igraph.Graph.Read_Edgelist(graph_file_list[i], directed = False)
        num_nodes = graph.vcount()

        for j in xrange(len(partition_file_lists[i])):
            
            partition_file = partition_file_lists[i][j]
            file_number_str = str(i + args.file_range_start) + '_' + str(j + 1)
            partition = partitionFromFile(partition_file, partition_separator)
            dotfile_path = args.output_path + name + '_' + dot_file_prefix + file_number_str + dot_file_suffix
            writeDotFile(graph_file_list[i], graph_separator, partition_file, partition_separator, dotfile_path)

            # To be written to CSV as program output
            result_lines = runIndividualGraphMetrics(graph, partition, dotfile_path, args.is_directed)

            print '\nSuccessfully analyzed individual metrics of ' + name + ' clustering for file number ' + file_number_str 

            writeRawDataLines(result_lines, name, str(num_nodes), str(i + args.file_range_start), str(j + 1))

            results.append(result_lines)

    return results

def goldStandardComparisonAnalysis(graph_file_list, partition_file_lists, gold_standard_file_lists, graph_separator, partition_separator, gold_standard_separator, name):
    """"""

    results = []
    
    for i in xrange(((len(graph_file_list) + len(partition_file_lists) + len(gold_standard_file_lists)) / 3)): # assume graph_file_list and partition_file_lists and gold_standard_file_lists are equal in length

        graph = igraph.Graph.Read_Edgelist(graph_file_list[i], directed = False)
        num_nodes = graph.vcount()

        for j in xrange(len(partition_file_lists[i])):

            partition_file = partition_file_lists[i][j]
            gold_standard_file = gold_standard_file_lists[i][0]
            file_number_str = str(i + args.file_range_start) + '_' + str(j + 1) 
            partition = partitionFromFile(partition_file, partition_separator)
            gold_standard = partitionFromFile(gold_standard_file, gold_standard_separator)

            result_lines = runComparisonGraphMetrics(graph, partition, gold_standard, args.is_directed)

            print '\nSuccessfully analyzed comparison metrics of ' + name + ' clustering for file number ' + file_number_str 
            
            writeRawDataLines(result_lines, name, str(num_nodes), str(i + args.file_range_start), str(j + 1))

            results.append(result_lines)

    return results

##############################
# Input Parameter Processing #
##############################

handleArgs()

raw_data_file_name = 'raw_data_s_' + str(args.file_range_start) + '_e_' + str(args.file_range_end) + '.csv'

graph_file_list = []
for i in range(args.file_range_start, args.file_range_end + 1):
    graph_file_list.append(args.graph_file_prefix + str(i) + args.graph_file_suffix)

if args.run_standard:
    gold_standard_file_list = []
    for i in range(args.file_range_start, args.file_range_end + 1):
        gold_standard_file_list.append([args.gold_standard_file_prefix + str(i) + args.gold_standard_file_suffix])

if len(args.clustering_file_suffixes) == 1:
    args.clustering_file_suffixes *= len(args.clustering_file_prefixes)

if len(args.clustering_file_separators) == 1:
    args.clustering_file_separators *= len(args.clustering_file_prefixes)

clustering_file_lists = []
for i in xrange(len(args.clustering_file_prefixes)):
    clustering_file_list = []
    prefix = args.clustering_file_prefixes[i]
    suffix = args.clustering_file_suffixes[i]
    for i in range(args.file_range_start, args.file_range_end + 1):
        same_file_clusterings = []
        for t in xrange(args.clusterings_per_graph):
            same_file_clusterings.append(prefix + str(i) + '_' + str(t+1) + suffix)
        clustering_file_list.append(same_file_clusterings)
    clustering_file_lists.append(clustering_file_list)

createPathIfNeeded(args.output_path)

##################
# Main execution #
##################

# Structural metrics on the graphs
#structure_lines = structuralMetricAnalysis(graph_file_list, args.graph_file_separator, 'network_structure')
#exportSummaryLines(structure_lines, 'network_structure')

# To collect clustering statistics
all_lines = []

# Individual metrics on "gold standard"
if args.run_standard:
    results = ['gold_standard']
    appendLines(individualMetricAnalysis(graph_file_list, gold_standard_file_list, args.graph_file_separator, args.gold_standard_file_separator, 'gold_standard'), results)
    all_lines.append(results)

# Individual and comparison metrics on input clusterings
for i in xrange(len(clustering_file_lists)):
    results = [args.clustering_file_names[i]]
    clustering_file_list = clustering_file_lists[i]
    clustering_file_separator = args.clustering_file_separators[i]
    appendLines(individualMetricAnalysis(graph_file_list, clustering_file_list, args.graph_file_separator, clustering_file_separator, args.clustering_file_names[i]), results)
    if args.run_standard:
        appendLines(goldStandardComparisonAnalysis(graph_file_list, clustering_file_list, gold_standard_file_list, args.graph_file_separator, clustering_file_separator, args.gold_standard_file_separator, args.clustering_file_names[i]), results)
    all_lines.append(results)

#exportSummaryStatistics(all_lines)

if 'oslom' in args.clustering_file_names:
    print '\nNote:\n\nOslom clustering produced overlapping partitions, but this code only considers the last assignment for each node in the clustering file.'

print('\nLogfile:')
for line in logfile_lines:
    print('\n' + line)
try:
    writeLinesToFile(logfile_lines, args.output_path + 'measure_logfile.log', mode = 'a')
except:
    print '\nError writing logfile_lines to file'
