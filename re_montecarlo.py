""""""

import numpy as np

H_residues = []
P_redidues = []
HP_residue = {
    **dict.fromkeys(H_residues, "H"),
    **dict.fromkeys(P_redidues, "P")
    }


def read_fasta(filename):
    """Read a fasta sequence and returns its sequence."""
    pass


def seq_to_hp(sequence):
    """Convert protein sequence to its corresponding HP sequence."""
    pass


if __name__ == "__main__":
    pass