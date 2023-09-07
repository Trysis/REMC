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
    """Returns positions adjacent to the chosen {position}."""
    moves = np.array([xy_up, xy_right, xy_left, xy_down], dtype=dtype)
    adjacent_pos = moves + np.array(position)
    return adjacent_pos.tolist()


def free_adjacent_positions(hp_coordinates, i, nonfree=False, dtype=np.int16):
    """Returns free (or non-free) position adjacent to residue {i}."""
    adjacent_pos = adjacent_positions(hp_coordinates[i], dtype=dtype)
    available_positions = []
    for adja_pos_i in adjacent_pos:
        if not nonfree and (adja_pos_i not in hp_coordinates):
            available_positions.append(adja_pos_i)
        if nonfree and (adja_pos_i in hp_coordinates):
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


def adjacent_neighbour(hp_coordinates, i, return_index=False):
    """Returns the coordinates (or index) of adjacent existant neighbours."""
    left = hp_coordinates[i-1].tolist() if i > 0 else None
    right = hp_coordinates[i+1].tolist() if i < (len(hp_coordinates) - 1) else None
    adjacent_pos = []
    if return_index:
        if left is not None:
            adjacent_pos.append(i-1)
        if right is not None:
            adjacent_pos.append(i+1)
    else:
        if left is not None:
            adjacent_pos.append(left)
        if right is not None:
            adjacent_pos.append(right)
    return adjacent_pos


def topological_neighbour(hp_coordinates, i):
    """Returns the coordinates of topological neighbours."""
    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    i = i - 1 if (i-1) in adjacent_nei_index else i
    # Mask to remove i-1 & i+1 neighbour
    mask = np.ones(len(hp_coordinates), dtype=bool)
    mask[np.array(adjacent_nei_index)] = 0
    # Coordinates without adjacent neighbour
    hp_coordinates = hp_coordinates[mask].tolist()
    topological_positions = free_adjacent_positions(hp_coordinates, i, nonfree=True)
    return topological_positions
    

def available_end_moves(hp_coordinates, i):
    """Returns available end moves positions from the first or last residue {i}."""
    if (i != 0) or (i != len(hp_coordinates)-1):
        return []
    adjacent_nei = adjacent_neighbour(hp_coordinates, i, return_index=True)
    adjacent_pos = free_adjacent_positions(hp_coordinates, adjacent_nei)
    return adjacent_pos


def available_pull_moves(hp_coordinates, i):
    pass


def plot_conformation(hp_coordinates):
    plt.scatter(hp_coordinates[:, 0], hp_coordinates[:, 1])
    plt.plot(hp_coordinates[:, 0], hp_coordinates[:, 1])
    min_xy = min([i[0] for i in hp_coordinates] + [i[1] for i in hp_coordinates])
    max_xy = max([i[0] for i in hp_coordinates] + [i[1] for i in hp_coordinates])

    #plt.xticks(np.arange(min_xy, max_xy, step=1))
    #plt.yticks(np.arange(min_xy, max_xy, step=1))
    plt.show()


if __name__ == "__main__":
    filename = "./data/A0A0C5B5G6.fasta"
    sequence = read_fasta(filename)
    hp_sequence = sequence_to_HP(sequence)
    hp_coordinates = initialize_coordinates(hp_sequence, random=True)
    for i in range(len(hp_coordinates)):
        print(topological_neighbour(hp_coordinates, i))
    plot_conformation(hp_coordinates)

