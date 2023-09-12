"""Replica Exchange Monte Carlo Algorithm (REMC)

This project was made for a course in the Master 2 Bioinformatics
course in the Université Paris Cité.

We needed to implement an algorithm from a chosen article pre-selected
by the teacher.

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

parser = argparse.ArgumentParser(
            prog = "REMC",
            description = "This program applies the Replica Exchange Monte " \
                          "Carlo Algorithm on a given conformation.",
            epilog = ""
            )

parser.add_argument('filepath', nargs='?', type=str, default="../data/A0A0C5B5G6.fasta", help="File path to fasta file.")
parser.add_argument('-o', '--image', type=str, default="../out/", help="Image output directory")
parser.add_argument('-n', '--replica', type=int, default=3, help="Replica number (please don't exagerate)")
parser.add_argument('-t', '--temperature', type=float, default=1.0, help="Temperature of the first replica")
parser.add_argument('-d', '--difference', type=float, default=0.0, help="Temperature difference from replica i to i+1")

args = parser.parse_args()
if not os.path.isfile(args.filepath):
    raise ValueError(f"Path is not correct : {args.filepath}")

if not os.path.isdir(args.image):
    raise ValueError(f"Path does not point to an existing directory : {args.image}")


if __name__ == "__main__":
    # TODO : Argparse
    filename = args.filepath    
    sequence = read_fasta(filename)
    print(sequence)
    exit()
    conf1 = Conformation(sequence, T=1, name="conf1")
    conf2 = Conformation(sequence, T=2, name="conf2")
    best_replica, conformations = \
            REMCSimulation(conformations=[conf1, conf2],
                           optimal_energy=-4, max_iter=10,
                           steps=50, neighbourhood_fct=vshd_neighbourhood,
                           move_on_step=True)

    best_replica.animate(save=True)
    print(f"{best_replica.E_list = }")
    print(f"{best_replica.T_list = }")
