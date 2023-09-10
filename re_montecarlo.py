""""""
import numpy as np
import matplotlib.pyplot as plt
import random

H_RES = ("V", "I", "F", "L", "M", "C", "W")
OTHER_RES = ("D", "E", "K", "R", "H", "Y", "S", "T", "N", "Q",
            "G", "A", "P",
            "B", "Z", "X", "J", "O", "U")  # specials

HP_RES = {
    **dict.fromkeys(H_RES, "H"),
    **dict.fromkeys(OTHER_RES, "P")
    }

cod_type = np.int16
# 2D
xy_up = np.array([0, 1], dtype=cod_type)
xy_right = np.array([1, 0], dtype=cod_type)
xy_left = np.array([-1, 0], dtype=cod_type)
xy_down = np.array([0, -1], dtype=cod_type)

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


def is_adjacent(position1, position2, dtype=np.int16):
    if position1[0] == position2[0]:
        return np.abs(position1[1] - position2[1]) == 1
    if position1[1] == position2[1]:
        return np.abs(position1[0] - position2[0]) == 1
    return False


def adjacent_positions(position, dtype=np.int16):
    """Returns positions adjacent to the chosen {position}."""
    moves = np.array([xy_up, xy_right, xy_left, xy_down], dtype=dtype)
    adjacent_pos = moves + position
    return adjacent_pos.tolist()

def diagonally_adjacent_positions(position, dtype=np.int16):
    """Returns diagonally adjacent position to the chosen {position}."""
    moves_left = np.array([xy_up, xy_down], dtype=dtype) + [-1, 0]
    moves_right = np.array([xy_up, xy_down], dtype=dtype) + [1, 0]
    diag_adjacent_pos = np.concatenate([moves_left, moves_right]) + position
    return diag_adjacent_pos

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

    # HP Coordinates without adjacent neighbour
    hp_coordinates_nonadja = hp_coordinates[mask].tolist()
    # Adjacent coordinates in list without i-1 & i+1
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
    # Free/Non-free selected positions
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
        return [], []
    
    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    available_adja_pos = free_adjacent_positions(hp_coordinates, adjacent_nei_index[0])
    return [i], available_adja_pos


def available_corner_moves(hp_coordinates, i):
    """Returns the available corner move position from a residue comprised in 1 to n-1."""
    if ((i == 0) or (i == len(hp_coordinates)-1)):
        return [], []

    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    available_adja_pos_left = free_adjacent_positions(hp_coordinates, adjacent_nei_index[0])
    available_adja_pos_right = free_adjacent_positions(hp_coordinates, adjacent_nei_index[1])
    # Search for free position mutually adjacent to i-1 & i+1 positions
    for available_pos in available_adja_pos_left:
        if available_pos in available_adja_pos_right:
            return [i], available_pos

    return [], []


def available_crank_shaft_moves(hp_coordinates, i):
    if ((i == 0) or (i == len(hp_coordinates)-1) or len(hp_coordinates)==3):
        return [], []

    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)  # i-1 & i+1
    adjacent_nei_left_idx = adjacent_neighbour(hp_coordinates, adjacent_nei_index[0], return_index=True)  # [i-2] & i
    adjacent_nei_right_idx = adjacent_neighbour(hp_coordinates, adjacent_nei_index[1], return_index=True)  # [i+2] & i

    adjacent_nei_left_idx.remove(i)  # i-2 or nothing
    adjacent_nei_right_idx.remove(i)  # i+2 or nothing
    
    # if either is True, then i is in a u-shaped conformation
    left, right = False, False
    to_move = []
    topological_nei_idx = []  # will contains either i-2 & i+1 or i-1 & i+2 or nothing
    if len(adjacent_nei_left_idx) == 1:  # i-2 exist
        topological_nei_left_idx = topological_neighbour(hp_coordinates, adjacent_nei_left_idx[0], return_index=True)
        for topological_idx in topological_nei_left_idx:
            # i-2 & i+1 are topological neighbours
            if topological_idx == adjacent_nei_index[1]:
                topological_nei_idx.append(adjacent_nei_left_idx[0])
                topological_nei_idx.append(adjacent_nei_index[1])
                to_move.append(i-1)
                to_move.append(i)
                left = True

    if len(adjacent_nei_right_idx) == 1:  # i+2 exist
        topological_nei_right_idx = topological_neighbour(hp_coordinates, adjacent_nei_right_idx[0], return_index=True)
        for topological_idx in topological_nei_right_idx:
            # i-1 & i+2 are topological neighbours
            if topological_idx == adjacent_nei_index[0]:
                topological_nei_idx.append(adjacent_nei_index[0])
                topological_nei_idx.append(adjacent_nei_right_idx[0])
                to_move.append(i)
                to_move.append(i+1)
                right = True

    if (left == False) & (right == False):
        return [], []
       
    free_adja_from_left = free_adjacent_positions(hp_coordinates, topological_nei_idx[0])
    free_adja_from_right = free_adjacent_positions(hp_coordinates, topological_nei_idx[1])
    
    # If an existing free position is available, the move can be performed
    for free_left in free_adja_from_left:
        for free_right in free_adja_from_right:
            if is_adjacent(free_left, free_right):
                return to_move, [free_left, free_right]
    
    return [], []
    

def available_pull_moves(hp_coordinates, i):
    """Returns the positions and residues to perform pull move on for residue {i}."""
    if (i == 0) or (i == len(hp_coordinates)-1):
        return [], []
    
    i_left = hp_coordinates[i-1].tolist()
    available_adja_pos = free_adjacent_positions(hp_coordinates, i) + [i_left]
    free_adja_from_right = free_adjacent_positions(hp_coordinates, i+1)
    print(f"{available_adja_pos = }")
    print(f"{free_adja_from_right = }")
    CL_positions = []
    for adja_i in available_adja_pos:
        for adja_right in free_adja_from_right:
            if is_adjacent(adja_i, adja_right):
                CL_positions.append([adja_i, adja_right])

    # Simplest case - Corner move - when i-1 in c position
    for c_position, _ in CL_positions:
        if i_left == c_position:
            return available_corner_moves(hp_coordinates, i)

    # No squares available (C & L positions), no moves can be performed
    if len(CL_positions) == 0:
        return [], []

    # In case if the movement can be performed on multiple direction
    selected_index = random_index(CL_positions)
    selected_CL = CL_positions[selected_index]
    available_pos, to_move = [selected_CL[1], selected_CL[0]], [i, i-1]
    j = i - 1
    while((j > 0) and not is_adjacent(available_pos[-1], hp_coordinates[j-1])):
        available_pos.append([hp_coordinates[j+1].tolist()])
        to_move.append(j-1)
        j = j - 1
    return to_move, available_pos

def is_H(converted_residue):
    return converted_residue == "H"

def conformation_energy(hp_sequence, hp_coordinates, dtype=np.int16, **kwargs):
    """"""
    if len(hp_sequence) != len(hp_coordinates):
        raise ValueError("arg1 and arg2 needs to have the same length.")
    if not isinstance(hp_coordinates, np.ndarray):
        hp_coordinates = np.array(hp_coordinates, dtype=dtype)
    # Penalty score
    H_penalty = kwargs.get("h_penalty", -1)
    otherwise_penalty = kwargs.get("o_penalty", 0)
    # 
    H_indices = [i for i, res in enumerate(hp_sequence) if is_H(res)]
    H_coordinates = hp_coordinates[H_indices, :].tolist()
    total_energy = 0
    n = len(hp_coordinates)
    for i, (h_i, h_coord_i) in enumerate(zip(H_indices[0: n-1], H_coordinates[0:n-1])):
        for h_j, h_coord_j in zip(H_indices[i+1:n], H_coordinates[i+1:n]):
            if np.abs(h_i - h_j) == 1:  # i&j are adjacent neighbour
                continue
            if is_adjacent(h_coord_i, h_coord_j):
                print(h_coord_i, h_coord_j)
                total_energy += H_penalty
            else:
                total_energy += otherwise_penalty

    print(f"{H_indices = }\n"
          f"{H_coordinates = }")
    return total_energy


def metropolis_criterion(hp_sequence_i, hp_coordinates_i,
                         hp_sequence_iplus1, hp_coordinates_iplus1, T,
                         **kwargs
):
    energy_i = conformation_energy(hp_sequence_i, hp_coordinates_i, **kwargs)
    energy_iplus1 = conformation_energy(hp_sequence_iplus1, hp_coordinates_iplus1, **kwargs)
    deltaEnergy = energy_iplus1 - energy_i
    if deltaEnergy <= 0:
        return 1.0
    else:
        np.exp((-deltaEnergy/T)).item()


def transition_probability():
    pass

def plot_conformation(hp_coordinates, hp_sequence=None):
    plt.scatter(hp_coordinates[:, 0], hp_coordinates[:, 1])
    if hp_sequence is not None:
        if(len(hp_coordinates) != len(hp_sequence)):
            raise ValueError("arg1 and arg2 needs to have the same length.")
        for i, hp_letter in enumerate(hp_sequence):
            plt.text(hp_coordinates[i, 0], hp_coordinates[i, 1], hp_letter)
    plt.plot(hp_coordinates[:, 0], hp_coordinates[:, 1])
    min_xy = min([i[0] for i in hp_coordinates] + [i[1] for i in hp_coordinates])
    max_xy = max([i[0] for i in hp_coordinates] + [i[1] for i in hp_coordinates])

    #plt.xticks(np.arange(min_xy, max_xy, step=1))
    #plt.yticks(np.arange(min_xy, max_xy, step=1))
    


if __name__ == "__main__":
    filename = "./data/A0A0C5B5G6.fasta"
    sequence = read_fasta(filename)
    hp_sequence = sequence_to_HP(sequence)
    hp_coordinates = initialize_coordinates(hp_sequence, random=True)
    print(conformation_energy(hp_sequence, hp_coordinates))
    plot_conformation(hp_coordinates, hp_sequence)
    plt.show()

