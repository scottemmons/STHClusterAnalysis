Installation
============

The brunt of installation work requires installing the package's dependencies (listed below). The package includes copies of all the dependencies that are not Python modules. All the dependencies besides the ModularityOptimizer and DEMON will need installation. Unzip all of the zip files, navigate to each resulting directory, and see the README's in each directory for installation instructions.

One additional step must be taken to configure the "binary_networks/" program. You must manually run the program once for it to produce the files necessary to be automatically run by this package's .py scripts. See the example installation workflow provided below for details.

Example Installation Workflow
=============================

1. Unzip binary_networks.tar.gz
2. Execute the terminal command `make` inside the binary_networks directory
3. Execute the terminal command `./benchmark -f flags.dat` followed by the terminal command `rm network.dat community.dat statistics.dat` inside the binary_networks directory in order manually to run the program once
4. Unzip clustering_programs.tar.gz
5. Execute the terminal command `./compile.sh` inside the clustering_programs_5_2 directory. Note that this will work from a Unix terminal, and that to run the program on Windows you can install MinGW from http://www.mingw.org/ (I haven't tried to install the program on Windows)
6. Unzip mutual3.tar.gz
7. Execute the terminal command `make` inside the mutual3 directory
8. Unzip gmap.zip
9. Execute the terminal command `make` inside the gmap/external/eba directory. Note that only this subset of the entire gmap program needs to be installed.
10. Additionally, you can test whether or not ModularityOptimizer.jar works with your installed version of Java by downloading the file karate_club_network.txt from http://www.ludowaltman.nl/slm/karate_club_network.txt, running the command `java -jar ModularityOptimizer.jar karate_club_network.txt karate_club_communities.txt 1 1.0 3 10 10 0 0`, and verifying that the output file karate_club_communities.txt in your directory matches the one found at http://www.ludowaltman.nl/slm/karate_club_communities.txt; see http://www.ludowaltman.nl/slm/ for more information about ModularityOptimizer.jar.

The above workflow details the specific steps required to install the external programs which this program calls. This program also requires Python 2.7 and the scientific programming libraries listed below in the Dependencies section. I recommend installing Anaconda (https://store.continuum.io/cshop/anaconda/), a Python distribution pre-configured with many, if not all, of the Python libraries this program uses.

Use
===

This network_cluster_analysis_package contains files which each correspond to a core function.
   - generate.py generates synthetic network graphs to be used as benchmarks for cluster analysis
   - cluster.py runs clustering algorithms over these graphs
   - measure.py clusters evaluation metrics over clusterings of a network graph
   - visualize.py visualizes the results of the cluster evaluation metrics

These files can be run individually via the command line for their respective functionalities.

The file "combined.py" combines the functionality of the individual files to run a complete cluster analysis workflow. For example, if you would like to run the full workflow for sizes of N = 1000 and N = 10000 and synthetic graphs of mixing parameter = 0.4, execute the command:
  python combined.py -n 1000 10000 --mu 0.4

Alternatively, you can specify a set of the individual methods to run with the "-m" flag. This assumes that the necessary files for each individual method exist in directories of the structure that "combined.py" would create them.
For example, you could run only the generate.py functionality on graphs of size N = 1000 and of mixing parameter = 0.4 by executing the command:
  python combined.py -n 1000 --mu 0.4 -m generate
Then, you could at a later time run only the cluster.py and measure.py functionalities on the graphs of size N = 1000 that were generated earlier by executing the command:
  python combined.py -m cluster measure -n 1000 --mu 0.4

The file incrementor.py facilitates executing instances of this program in parallel by incrementing the seed for the random number generator in the binary_networks program. It should be called after executing the generate functionality of this code to increment the seed for later executions. If files are named and located according to installation defaults, no arguments to this script should be required, and it can be run with the command:
  python incrementor.py
  
Additionally, the start and end flags can be used to specify which trial numbers to execute. The defaults are start = 1 and end = 10, but the entire workflow could be run for only one trial by calling the command:
  python combined.py -n 1000 --mu 0.4
This functionality can be used to run multiple instances of the experiment in parallel. For example, to run three trials of the code in parallel, execute in parallel these three commands:
  python combined.py -n 1000 --mu 0.4 -s 1 -e 1
  python combined.py -n 1000 --mu 0.4 -s 2 -e 2
  python combined.py -n 1000 --mu 0.4 -s 3 -e 3
  
Example Use Workflow
====================

Here is a workflow for network graphs of size N = 1000 and mixing parameter = 0.4 broken down into its individual components:
1. Execute 'python combined.py -n 1000 --mu 0.4 -m generate'. This will generate synthetic graphs in the folder generated_benches/n_1000
2. Execute 'python combined.py -n 1000 --mu 0.4 -m cluster'. This will cluster the generated synthetic graphs, producing clustering files in the folder generated_benches/n_1000
3. Execute 'python combined.py -n 1000 --mu 0.4 -m measure'. This will measure the properties, using both gold standard comparison and stand-alone metrics, of the clusterings of the network graphs, producing measurement values of the prefix 'raw_data' in the folder generated_benches/n_1000
4. Execute 'python combined.py -n 1000 --mu 0.4 -m visualize'. This will create visualizations in the form of line graphs and box-and-whisker plots of the measurements produced by measure.py. The visualizations can be found in the generated_visualizations folder
5. Execute 'python incrementor.py' to increment the random number seed of the binary_networks program for future trials.

Functionality identical to the previous workflow can be achieved by using the -m flag to combine calls to combined.py:
1. Execute 'python combined.py -n 1000 --mu 0.4 -m generate cluster measure visualize'. This achieves steps 1-4 from the previous workflow.
2. Execute 'python incrementor.py' to increment the random number seed of the binary_networks program for future trials.
Note that as long as generate, cluster, measure, and visualize functionalities are called in that order, the program will execute correctly. Consequently, there are a variety of ways in which combined.py can be called to run these experiments, giving you control over which aspects of the experiments are run at any given time.

You could run fifty trials of the experiments for graphs of sizes N = 1000, N = 10000, and N = 100000 and mixing parameter = 0.4 with these commands:
1. python combined.py -n 1000 10000 100000 --mu 0.4 -m generate cluster measure visualize -s 1 -e 50
2. python incrementor.py
Note that visualize.py produces visualizations showing both how the values of one metric change with the size of the graph and comparing different metric values on one graph size.

Dependencies
============

The entire package is built on Python 2.7, which is required. In addition, individual files have the following dependencies:

generate.py
   - Lancich. benchmark generation "binary_networks/" @ https://sites.google.com/site/andrealancichinetti/files/binary_networks.tar.gz

cluster.py
   - Lancich. clustering repository "clustering_programs_5_2/" @ https://sites.google.com/site/andrealancichinetti/clustering_programs.tar.gz (Note that the version of this program in the .tar.gz file I shared is a slightly modified version of the original. Specifically, I modified the select.py and methods.py files in order not to run the 'cvis' step, which is auxiliary to the clusteirng functionality of the code and adds to its runtime.)
   - Leiden clustering jar "ModularityOptimization.jar" @ http://www.ludowaltman.nl/slm/ModularityOptimizer.jar
      - Java (i.e. JRE)
   - DEMON clustering archive, which must be unzipped into "demon_py/" @ http://www.michelecoscia.com/wp-content/uploads/2013/07/demon_py.zip
      - python module networkx

measure.py
   - python module igraph
   - python module sklearn
   - my fork of GMap "~/Documents/gmap/external/eba" @ https://github.com/scottemmons/gmap.git
   - Lancich. NMI "mutual3/" @ https://sites.google.com/site/andrealancichinetti/mutual3.tar.gz

visualize.py
   - python module pandas
   - python module numpy
   - python module matplotlib (specifically matplotlib.pyplot)

Notes
=====

Right now the combined.py code is set up to handle the random number seed correctly for 100 trials at the sizes N = 1000, 10000, 100000, and 1000000. If you want to run this experiment in a different way, look over how the random number seed is being handled and update the code accordingly to facilitate your workflow.

"Blondel" appears in the code as another name for the Louvain algorithm, but it does not appear in the paper in order not to imply undue credit to Blondel over his collaborators.
