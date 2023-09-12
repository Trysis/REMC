"""Replica Exchange Monte Carlo Algorithm (REMC)

This project was made for a course in the Master 2 Bioinformatics
course in the Université Paris Cité.

We needed to implement an algorithm from a chosen article pre-selected
by the teacher.

Based on :
[Thachuk C, Shmygelska A, Hoos HH.
A replica exchange Monte Carlo algorithm for protein folding in the HP model. BMC Bioinformatics.
2007 Sep 17;8:342. doi: 10.1186/1471-2105-8-342. PMID: 17875212; PMCID: PMC2071922.]

"""

import os

# Argument parser
import argparse

# Local module
from conformation import Conformation
from utility import *
from moves import *

__authors__ = "Roude JEAN MARIE"
__contact__ = ("roude.etu@gmail.com", "roude.bioinfo@gmail.com")
__date__ = "01-09-2023"
__version__ = "1.0.0"

neighbourhood_algorithm = {"vshd": vshd_neighbourhood,
                           "pull": pull_move_neighbourhood,
                           "hybrid": hybrid_neighbourhood}

parser = argparse.ArgumentParser(
            prog = "REMC",
            description = "This program applies the Replica Exchange Monte " \
                          "Carlo Algorithm on a given conformation.",
            epilog = ""
            )

parser.add_argument('filepath', nargs='?', type=str, default="../data/A0A0C5B5G6.fasta", help="File path to fasta file.")
parser.add_argument('-o', '--output', type=str, default="../out/", help="output directory")
parser.add_argument('-n', '--replica', type=int, default=3, help="Replica number (please don't exagerate)")
parser.add_argument('-b', '--random', action='store_true', default=False, help="Should the conformation be linearly initialized or randomly ?")
# Temperature
parser.add_argument('-t', '--temperature', type=float, default=1.0, help="Temperature of the first replica")
parser.add_argument('-d', '--difference', type=float, default=0.0, help="Temperature difference from replica i to i+1")
# Algorithm parameters
parser.add_argument('-a', '--algorithm', type=str.lower, choices=['vshd', 'pull', 'hybrid'], default="vshd", help="Chosen neigbourhood algorithm between vshd, pull and hybrid (both)")
parser.add_argument('-s', '--steps', type=int, default=20, help="Number of steps for the Monte Carlo search algorithm")
parser.add_argument('-r', '--rsteps', type=int, default=10, help="Number of times the Replica Exchange Algorithm will be performed")
parser.add_argument('-e', '--energy', type=float, default=-99999, help="Optimal energy we want our conformation to stop on (if steps*rsteps not reached)")
parser.add_argument('-p', '--proba', type=float, default=0.5, help="Hybrid probability algorithm when chosen.")

# Get args
args = parser.parse_args()

# Exceptions
if not os.path.isfile(args.filepath):
    raise ValueError(f"Path is not correct : {args.filepath}")

if not os.path.isdir(args.output):
    raise ValueError(f"Path does not point to an existing directory : {args.image}")

if args.replica <= 0:
    raise ValueError(f"Number of replica should be strictly positive (>0). option -n")

if args.temperature <= 0:
    raise ValueError(f"Temperature should be strictly positive. option -t")

if args.difference < 0:
    raise ValueError(f"Temperature differences between replicas should be positive. option -d")

if args.steps <= 0:
    raise ValueError(f"Number steps should be strictly positive. option -s")

if args.rsteps <= 0:
    raise ValueError(f"Number of REMC steps should be strictly positive. option -r")

if args.algorithm == "hybrid":
    neighbourhood_algorithm[args.algorithm] = neighbourhood_algorithm[args.algorithm](args.proba)

if (args.proba < 0) or (args.proba > 1):
    raise ValueError(f"Probability should be comprised between 0 and 1. option -p")

# Retrieve args
output_directory = args.output
# Conformation(s)
n = args.replica
random_conformation = args.random
# Temperature
temperature = args.temperature
temperature_differences = args.difference
temperature_set = [temperature+temperature_differences*i for i in range(n)]
# Algorithm
nei_algorithm = neighbourhood_algorithm[args.algorithm]
steps = args.steps
rsteps = args.rsteps
optimal_e = args.energy

if __name__ == "__main__":
    filename = args.filepath    
    sequence = read_fasta(filename)
    conformations = [Conformation(sequence, T=temp, name=f"r{i}_len={len(sequence)}", random=random_conformation) \
                     for i, temp in enumerate(temperature_set)]

    best_replica, conformations = \
            REMCSimulation(conformations=conformations,
                           optimal_energy=optimal_e, max_iter=rsteps,
                           steps=steps, neighbourhood_fct=nei_algorithm,
                           move_on_step=True)

    for replica in conformations:
        replica.animate(show=False, save=True, saveto=output_directory)

    best_replica.animate(show=True, save=True, saveto=output_directory, filename="best_replica")
    best_replica.plot(initial=True, show=False, save=True, saveto=output_directory, filename="best_r_initial")
    best_replica.plot(show=False, save=True, saveto=output_directory, filename="best_r")
    best_replica.plot_energy(save=True, saveto=output_directory)