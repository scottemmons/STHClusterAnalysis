Installation
============

The brunt of installation work requires installing the package's dependencies (listed below). The package includes copies of all the dependencies that are not Python modules. All the dependencies besides the ModularityOptimizer will need installation. Unzip all of the zip files, navigate to each resulting directory, and see the README's in each directory for installation instructions.

Use
===

This network_cluster_analysis_package contains files which each correspond to a core function.
   - generate.py will generate synthetic network graphs to be used as benchmarks for cluster analysis
   - cluster.py will run clustering algorithms over these graphs
   - measure.py run cluster evaluation metrics over clusterings of a network graph
   - visualize.py visualizes the results of the cluster evaluation metrics

These files can be run individually via the command line for their respective functionalities.

The file "combined.py" combines the functionality of the individual files to run a complete cluster analysis workflow. For example, if you would like to run the full workflow for sizes of N = 100, N = 1000, and N = 10000, execute the command:
  python combined.py -n 100 1000 10000
Alternatively, you can specify a set of the individual methods to run with the "-m" flag. This assumes that the necessary files for each individual method exist in directories of the structure that "combined.py" would create them.
For example, you could run only the generate.py functionality on graphs of size N = 100 and N = 1000 by executing the command:
  python combined.py -n 100 1000 -m generate
Then, you could at a later time run only the cluster.py and measure.py functionalities on the graphs of size N = 100 and N = 1000 that were generated earlier by executing the command:
  python combined.py -m cluster measure -n 100 1000

Dependencies
============

The entire package is built on Python 2.7, which is required. In addition, individual files have the following dependencies:

generate.py
   - Lancich. benchmark generation "binary_networks/" @ https://sites.google.com/site/andrealancichinetti/files/binary_networks.tar.gz

cluster.py
   - Lancich. clustering repository "clustering_programs_5_2/" @ https://sites.google.com/site/andrealancichinetti/clustering_programs.tar.gz
   - Leiden clustering jar "ModularityOptimization.jar" @ http://www.ludowaltman.nl/slm/ModularityOptimizer.jar
      - Java (i.e. JRE)

measure.py
   - python module igraph
   - python module sklearn
   - my fork of GMap "~/Documents/gmap/external/eba" @ https://github.com/scottemmons/gmap.git
   - Lancich. NMI "mutual3/" @ https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz

visualize.py
   - python module pandas
   - python module numpy
   - python module matplotlib (specifically matplotlib.pyplot)
