import argparse
import os
import shutil

package_name = "network_cluster_analysis_package"
package_path = os.path.join(os.getcwd(), package_name)
package_components = ['combined.py', 'generate.py', 'cluster.py', 'measure.py', 'visualize.py', 'incrementor.py']

def handleArgs():
    """Handles the command-line input arguments, placing them in the global Namespace variable 'args'."""
    parser = argparse.ArgumentParser(description="Automates the command line calls necessary to execute the full clustering analysis script workflow")
    source_destination_group = parser.add_mutually_exclusive_group(required=True)
    source_destination_group.add_argument("--to", action="store_true", help="indicates that files should be moved to the package source", dest="move_to")
    source_destination_group.add_argument("--from", action="store_true", help="indicates that files should be moved from the package source", dest = "move_from")
    
    global args
    args = parser.parse_args()

handleArgs()

for component in package_components:
    cwd = os.path.join(os.getcwd(), component)
    pkg = os.path.join(package_path, component)
    
    if args.move_to:
        shutil.copyfile(cwd, pkg)
    elif args.move_from:
        shutil.copyfile(pkg, cwd)
