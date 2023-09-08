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


def initialize_coordinates(hp_sequence, random=False, dtype=np.int16):
    """Returns a set of coordinates for the associated sequence."""
    positions = [[0, 0]]  # position of first residue
    # Coordinates are aligned
    if not random:
        for i in range(1, len(hp_sequence)):
            positions.append([positions[0][0] + i, positions[0][1]])
    # The path chosen is to be set at random
    elif random:
        i = 0
        moves = [[]] * (len(hp_sequence) - 1)
        # While residues are not all placed
        while(len(positions) != len(hp_sequence)):
            available_positions = free_adjacent_positions(positions, i, dtype=dtype)
            # If path is blocked we backtrack
            if len(available_positions) == 0:
                while(len(moves[i-1]) == 0):
                    i = i - 1
                available_positions = moves[i-1]
                positions = positions[:i]
                i = i - 1
            # Else we define the residue position
            selected_index = random_index(available_positions)
            selected_position = available_positions.pop(selected_index)
            positions.append(selected_position)
            moves[i] = available_positions
            i = i + 1

    return np.array(positions, dtype=dtype)


def random_index(positions):
    """Return a random index from the available indices."""
    selected_index = random.randint(0, len(positions) - 1)
    return selected_index


def rotate_position(position, clockwise=True):
    """Perform a 90 degree rotation to the specified position."""
    x, y = position[0], position[1]
    if clockwise:
        position[0], position[1] = y, -x
    elif not clockwise:
        position[0], position[1] = -y, x
    return position


def adjacent_positions(position, dtype=np.int16):
    """Returns positions adjacent to the chosen {position}."""
    moves = np.array([xy_up, xy_right, xy_left, xy_down], dtype=dtype)
    adjacent_pos = moves + np.array(position)
    return adjacent_pos.tolist()


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


def topological_neighbour(hp_coordinates, i, return_index=False):
    """Returns the coordinates of topological neighbours."""
    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    i_shifted = i - 1 if (i-1) in adjacent_nei_index else i

    # Mask to remove i-1 & i+1 neighbour
    mask = np.ones(len(hp_coordinates), dtype=bool)  # full of ones
    mask[np.array(adjacent_nei_index)] = 0

    # Coordinates without adjacent neighbour
    hp_coordinates_nonadja = hp_coordinates[mask].tolist()
    topological_positions = free_adjacent_positions(hp_coordinates_nonadja, i_shifted, nonfree=True)

    if return_index:
        indices = []
        for idx, hp_position in enumerate(hp_coordinates.tolist()):
            if hp_position in topological_positions:
                indices.append(idx)
        return indices

    return topological_positions


def free_adjacent_positions(hp_coordinates, i, nonfree=False, dtype=np.int16):
    """Returns free (or non-free) position adjacent to residue {i}."""
    if isinstance(hp_coordinates, np.ndarray):
        hp_coordinates = hp_coordinates.tolist()
    if not isinstance(hp_coordinates, list):
        hp_coordinates = list(hp_coordinates)

    adjacent_pos = adjacent_positions(hp_coordinates[i], dtype=dtype)
    # 
    selected_positions = []
    for adja_pos_i in adjacent_pos:
        if (not nonfree) and (adja_pos_i not in hp_coordinates):
            selected_positions.append(adja_pos_i)
        if nonfree and (adja_pos_i in hp_coordinates):
            selected_positions.append(adja_pos_i)

    return selected_positions
    

def available_end_moves(hp_coordinates, i):
    """Returns available end moves positions from the first or last residue {i}."""
    if not ((i == 0) or (i == len(hp_coordinates)-1)):
        return []
    
    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    adjacent_pos = free_adjacent_positions(hp_coordinates, adjacent_nei_index[0])
    return adjacent_pos


def available_pull_moves(hp_coordinates, i):
    """Returns the available pull move position from a residue comprised in 1 to n-1."""
    if ((i == 0) or (i == len(hp_coordinates)-1)):
        return []

    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    available_pos_left = free_adjacent_positions(hp_coordinates, adjacent_nei_index[0])
    available_pos_right = free_adjacent_positions(hp_coordinates, adjacent_nei_index[1])
    # Search for free position mutually adjacent to i-1 & i+1 positions
    for available_pos in available_pos_left:
        if available_pos in available_pos_right:
            return available_pos

    return []


def available_crank_shaft_moves(hp_coordinates, i):
    if ((i == 0) or (i == len(hp_coordinates)-1) or len(hp_coordinates)==3):
        return []

    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)  # i-1 & i+1
    adjacent_nei_left_idx = adjacent_neighbour(hp_coordinates, adjacent_nei_index[0], return_index=True)  # [i-2] & i
    adjacent_nei_right_idx = adjacent_neighbour(hp_coordinates, adjacent_nei_index[1], return_index=True)  # [i+2] & i

    adjacent_nei_left_idx.remove(i)  # i-2 or nothing
    adjacent_nei_right_idx.remove(i)  # i+2 or nothing

    # Case if both i-1 & i+1 don't have neighbours
    if len(adjacent_nei_left_idx + adjacent_nei_right_idx) == 0:
        return []
    
    u_shaped = False
    topological_nei_idx = []
    if len(adjacent_nei_left_idx) == 1:  # i-2 exist
        topological_nei_left_idx = topological_neighbour(hp_coordinates, adjacent_nei_left_idx[0], return_index=True)
        for topological_idx in topological_nei_left_idx:
            # i-2 & i+1 are topological neighbours
            if topological_idx == adjacent_nei_index[1]:
                topological_nei_idx.append(adjacent_nei_left_idx[0], adjacent_nei_index[1])
                u_shaped = True
                break

    elif len(adjacent_nei_right_idx) == 1:  # i+2 exist
        topological_nei_right_idx = topological_neighbour(hp_coordinates, adjacent_nei_right_idx[0], return_index=True)
        for topological_idx in topological_nei_right_idx:
            # i-1 & i+2 are topological neighbours
            if topological_idx == adjacent_nei_index[0]:
                topological_nei_idx.append(adjacent_nei_index[0], adjacent_nei_right_idx[0])
                u_shaped = True
                break
    
    if not u_shaped:
        return []
    
    # free positions
    
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
        print(available_pull_moves(hp_coordinates, i))
    print(available_end_moves(hp_coordinates, len(hp_coordinates)-1))
    print()
    print(topological_neighbour(hp_coordinates, len(hp_coordinates)-1, return_index=False))
    print(topological_neighbour(hp_coordinates, len(hp_coordinates)-1, return_index=True))
    plot_conformation(hp_coordinates)

