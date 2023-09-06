""""""
import numpy as np
import matplotlib.pyplot as plt
import random

H_RES = ["V", "I", "F", "L", "M", "C", "W"]
OTHER_RES = ["D", "E", "K", "R", "H", "Y", "S", "T", "N", "Q"] + \
            ["G", "A", "P"] + \
            ["B", "Z", "X", "J", "O", "U"]  # specials

HP_RES = {
    **dict.fromkeys(H_RES, "H"),
    **dict.fromkeys(OTHER_RES, "P")
    }

# 2D
xy_up = [0, 1]
xy_right = [1, 0]
xy_left = [-1, 0]
xy_down = [0, -1]


def read_fasta(filename):
    """Read the first fasta sequence and returns its sequence."""
    sequence = None
    with open(filename, "r") as fasta_in:
        header = fasta_in.readline().strip()
        sequence = fasta_in.read().strip()
        last_index = sequence.find(">")
        if last_index != -1:
            sequence = sequence[:last_index]
        if header[0] != ">":
            sequence = header + sequence

    return sequence


def sequence_to_HP(sequence):
    """Convert protein sequence to its corresponding HP sequence."""
    sequence = sequence.upper()
    hp_sequence = ""
    for res in sequence:
        hp_sequence += HP_RES[res]
    return hp_sequence


def random_index(positions):
    """Return a random index from the available indices."""
    selected_index = random.randint(0, len(positions) - 1)
    return selected_index


def adjacent_positions(position, dtype=np.int16):
    moves = np.array([xy_up, xy_right, xy_left, xy_down], dtype=dtype)
    adjacent_positions = moves + np.array(position)
    return adjacent_positions.tolist()

def free_adjacent_positions(hp_coordinates, i, dtype=np.int16):
    """Returns free position adjacent to residue i."""
    adjacent_positions = adjacent_positions(hp_coordinates[i], dtype=dtype)
    available_positions = []
    for adja_pos_i in adjacent_positions:
        if adja_pos_i not in hp_coordinates.tolist():
            available_positions.append(adja_pos_i)

    return available_positions


def topological_neighbour(hp_coordinates, i, dtype=np.int16):
    pass

def initialize_coordinates(hp_sequence, random=False, dtype=np.int16):
    positions = [[0, 0]]  # position of first residue
    if not random:
        for i in range(1, len(hp_sequence)):
            positions.append([i, 0])
    else:  # If path is blocked before adding all residue, it will fail
        # TODO : Parcours en largeur, ou profondeur avec choix random (selon la taille de la séquence)
        for i in range(1, len(hp_sequence)):
            np_positions = np.array(positions, dtype=dtype)
            available_position = free_adjacent_positions(np_positions, i-1, dtype=dtype)
            select_index = random_index(available_position)
            selected_position = available_position[select_index]
            positions.append(selected_position)

    return np.array(positions, dtype=dtype)


def plot_conformation(hp_coordinates):
    plt.scatter(hp_coordinates[:, 0], hp_coordinates[:, 1])
    plt.plot(hp_coordinates[:, 0], hp_coordinates[:, 1])
    plt.show()


if __name__ == "__main__":
    filename = "./data/A0A0C5B5G6.fasta"
    sequence = read_fasta(filename)
    hp_sequence = sequence_to_HP(sequence)
    hp_coordinates = initialize_coordinates(hp_sequence, random=False)
    plot_conformation(hp_coordinates)

