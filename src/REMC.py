"""Replica Exchange Monte Carlo Algorithm (REMC)

This project was made for a course in the Master 2 Bioinformatics
in the Université Paris Cité.

We needed to implement an algorithm from a chosen article pre-selected
by the teacher.

"""

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

if __name__ == "__main__":
    # TODO : Argparse
    filename = "../data/A0A0C5B5G6.fasta"
    sequence = read_fasta(filename)
    conf1 = Conformation(sequence, T=1, name="conf1")
    conf2 = Conformation(sequence, T=1, name="conf2")
    best_replica, conformations = \
            REMCSimulation(conformations=[conf1, conf2],
                           optimal_energy=-4, max_iter=10,
                           steps=50, neighbourhood_fct=vshd_neighbourhood,
                           move_on_step=True)

    best_replica.animate(save=True)