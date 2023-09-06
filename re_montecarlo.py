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
    adjacent_pos = adjacent_positions(hp_coordinates[i], dtype=dtype)
    available_positions = []
    for adja_pos_i in adjacent_pos:
        if adja_pos_i not in hp_coordinates:
            available_positions.append(adja_pos_i)

    return available_positions


def initialize_coordinates(hp_sequence, random=False, dtype=np.int16):
    positions = [[0, 0]]  # position of first residue
    # Coordinates are aligned
    if not random:
        for i in range(1, len(hp_sequence)):
            positions.append([positions[0][0] + i, positions[0][1]])
    # The path chosen is to be set at random
    elif random:
        i = 0
        moves = [[]] * (len(hp_sequence) - 1)
        while(len(positions) != len(hp_sequence)):
            available_positions = free_adjacent_positions(positions, i, dtype=dtype)
            # If path is blocked we backtrack
            if len(available_positions) == 0:
                while(len(moves[i-1]) == 0):
                    i = i - 1
                available_positions = moves[i-1]
                positions = positions[:i]
                i = i - 1
            # Else path is selected
            selected_index = random_index(available_positions)
            selected_position = available_positions.pop(selected_index)
            positions.append(selected_position)
            moves[i] = available_positions
            i = i + 1

    return np.array(positions, dtype=dtype)


def adjacent_neighbour(hp_coordinates, i, dtype=np.int16):
    left = hp_coordinates[i-1] if i > 0 else None
    right = hp_coordinates[i+1] if i < (len(hp_coordinates) - 1) else None
    if left is None:
        return [right]
    if right is None:
        return [left]
    return [left, right]

def topological_neighbour(hp_coordinates, i, dtype=np.int16):
    


def plot_conformation(hp_coordinates):
    plt.scatter(hp_coordinates[:, 0], hp_coordinates[:, 1])
    plt.plot(hp_coordinates[:, 0], hp_coordinates[:, 1])
    plt.show()


if __name__ == "__main__":
    filename = "./data/A0A0C5B5G6.fasta"
    sequence = read_fasta(filename)
    hp_sequence = sequence_to_HP(sequence)
    hp_coordinates = initialize_coordinates(hp_sequence, random=True)
    plot_conformation(hp_coordinates)
    print()

