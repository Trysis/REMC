""""""

# Python default modules
import random

# Data gestion
import numpy as np
import matplotlib.pyplot as plt

# Coordinates type
COORD_TYPE = np.int16

# 2D
XY_UP = np.array([0, 1], dtype=COORD_TYPE)
XY_RIGHT = np.array([1, 0], dtype=COORD_TYPE)
XY_LEFT = np.array([-1, 0], dtype=COORD_TYPE)
XY_DOWN = np.array([0, -1], dtype=COORD_TYPE)


def is_adjacent(position1, position2):
    """Check if two 2D positions are adjacent.

    position1, position2: list or numpy.ndarray of shape (2,)
        array-like containing x and y coordinates at i=0 & i=1
        indices.

    Returns: bool
        boolean indicating if the two positions are adjacent
        or not.

    """
    if position1[0] == position2[0]:
        return np.abs(position1[1] - position2[1]) == 1
    if position1[1] == position2[1]:
        return np.abs(position1[0] - position2[0]) == 1
    return False


def adjacent_positions(position, dtype=np.int16):
    """Returns adjacent positions to the chosen {position}.

    position: list or (numpy.ndarray -> of shape (2,))
        array-like containing x and y coordinates at i=0 & i=1
        indices.

    dtype: type
        numpy array type for the coordinates

    Returns: list -> of shape (4, 2)
        Array containing at each index the adjacent position to
        {position}.

    """
    moves = np.array([XY_UP, XY_RIGHT, XY_LEFT, XY_DOWN], dtype=dtype)
    adjacent_pos = moves + position
    return adjacent_pos.tolist()


def adjacent_neighbour(hp_coordinates, i, return_index=False):
    """Returns the associated coordinates or index adjacent to residue {i}.
    
    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get adjacent neighbour on.

    return_index: bool
        If True the functions return the indices of the adjacent neighbour
        residue.
    
    Returns: list -> of shape (x, 2) for residue coordinates or (x ,1) for indices
        Array containing the coordinates or indices of the adjacent neighbour of
        residue {i}. Adjacent neighbour correspond to existing {i-1} & {i+1} residue.

    """
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
    """Returns the coordinates of topological neighbours.

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get topological neighbour on.

    return_index: bool
        If True the functions return the indices of the adjacent neighbour
        residue.

    Returns: list -> of shape (x, 2) for residue coordinates or (x ,1) for indices
        Array containing the coordinates or indices of the topological neighbour of
        residue {i}. Topological neighbour correspond to residue that are not adjacent
        neighbour, but have coordinates adjacent to the residue {i}.

    """
    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    i_shifted = i - 1 if (i-1) in adjacent_nei_index else i

    # Mask to remove i-1 & i+1 neighbour
    mask = np.ones(len(hp_coordinates), dtype=bool)  # full of True
    mask[np.array(adjacent_nei_index)] = 0  # False for i-1 & i+1 index

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
    """Returns free (or non-free) position adjacent to residue {i}.

    hp_coordinates: array-like or numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get free/occupied positions on.

    nonfree: bool
        When true the function will return the occupied adjacent coordinates,
        otherwise it will returns the free adjacent coordinates to residue {i}.

    dtype: type
        numpy array type for the coordinates

    Returns: list -> of shape (x, 2)
        Array containing the coordinates of the available (or occupied) adjacent
        positions to the residue {i}.

    """
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
    """Returns an available end move position.
        Can only be performed on the first or last residue {i} in (0, n-1).

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    Returns: list, list -> of shape: (x, 1), (x, 2)
        Returns two array, the first contains the residue that is
        to be moved and the second contains the coordinates where residue
        in the first array can move.

    """
    if not ((i == 0) or (i == len(hp_coordinates)-1)):
        return [], []
    
    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    available_adja_pos = free_adjacent_positions(hp_coordinates, adjacent_nei_index[0])
    # If no existing end move position are available, then return empty lists
    if len(available_adja_pos) == 0:
        return [], []

    return [i], [random.choice(available_adja_pos)]


def available_corner_moves(hp_coordinates, i):
    """Returns an available corner move position.
        Residue {i} is comprised in [1; n-1].

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    Returns: list, list -> of shape: (x, 1), (x, 2)
        Returns two array, the first contains the residue that is
        to be moved and the second contains the coordinates where residue
        in the first array can move.

    """
    if ((i == 0) or (i == len(hp_coordinates)-1)):
        return [], []

    adjacent_nei_index = adjacent_neighbour(hp_coordinates, i, return_index=True)
    available_adja_pos_left = free_adjacent_positions(hp_coordinates, adjacent_nei_index[0])
    available_adja_pos_right = free_adjacent_positions(hp_coordinates, adjacent_nei_index[1])
    # Search for free position mutually adjacent to i-1 & i+1 positions
    for available_pos in available_adja_pos_left:
        if available_pos in available_adja_pos_right:
            return [i], [available_pos]

    return [], []


def available_crank_shaft_moves(hp_coordinates, i):
    """Returns an available crank shaft move position.
        {i} should be comprised in [1, n-1] and
        A sequence of length 3 or less cannot perform this move.
    
    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    Returns: list, list -> of shape: (x, 1), (x, 2)
        Returns two array, the first contains the residue that is
        to be moved and the second contains the coordinates where residue
        in the first array can move.

    """
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
    while((j > 0) and not is_adjacent(available_pos[-1], hp_coordinates[j-1].tolist())):
        available_pos.append(hp_coordinates[j+1].tolist())
        to_move.append(j-1)
        j = j - 1
    return to_move, available_pos


def vshd_neighbourhood(hp_coordinates, i):
    available_movements = {
        "end": available_end_moves(hp_coordinates, i),
        "corner": available_corner_moves(hp_coordinates, i),
        "crank_shaft": available_crank_shaft_moves(hp_coordinates, i)
    }
    # Delete non-performable moves
    to_del = []
    for key, (res_indices, res_positions) in available_movements.items():
        if len(res_indices) == 0:  # No move available for the specific key
            to_del.append(key)

    for key in to_del:
        del available_movements[key]

    # Returns available moves
    return available_movements


def pull_move_neighbourhood(hp_coordinates, i):
    available_movements = {
        "pull_move": available_pull_moves(hp_coordinates, i),
    }
    # Delete non-performable moves
    to_del = []
    for key, (res_indices, res_positions) in available_movements.items():
        if len(res_indices) == 0:  # No move available for the specific key
            to_del.append(key)

    for key in to_del:
        del available_movements[key]

    # Returns available moves
    return available_movements


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
                total_energy += H_penalty
            else:
                total_energy += otherwise_penalty

    return total_energy


def metropolis_criterion(hp_sequence, hp_coordinates_i, hp_coordinates_iplus1, T, **kwargs):
    """"""
    energy_i = conformation_energy(hp_sequence, hp_coordinates_i, **kwargs)
    energy_iplus1 = conformation_energy(hp_sequence, hp_coordinates_iplus1, **kwargs)
    deltaEnergy = energy_iplus1 - energy_i
    if deltaEnergy <= 0:
        return 1.0
    else:
        return np.exp((-deltaEnergy/T)).item()


def re_criterion(hp_sequence, hp_coordinates_i, T_i, hp_coordinates_j, T_j, **kwargs):
    """Replica exchange probability between two conformations."""
    energy_i = conformation_energy(hp_sequence, hp_coordinates_i, **kwargs)
    energy_j = conformation_energy(hp_sequence, hp_coordinates_j, **kwargs)
    criterion = ((1/T_j) - (1/T_i)) * (energy_i - energy_j)
    if criterion <= 0:
        return 1.0
    else:
        return np.exp(-criterion).item()


def mc_move(hp_coordinates, i, neighbourhood_fct=vshd_neighbourhood):
    """"""
    available_movements = neighbourhood_fct(hp_coordinates, i)
    if len(available_movements) == 0:
        return False, hp_coordinates

    # Change coordinates based on movement
    hp_coordinates_iplus1 = np.copy(hp_coordinates)
    move = random.choice(list(available_movements.keys()))
    for res_index, res_new_pos in zip(available_movements[move][0], available_movements[move][1]):
        hp_coordinates_iplus1[res_index] = res_new_pos

    return True, hp_coordinates_iplus1


def MCsearch(hp_sequence, hp_coordinates, T,
             steps=100, neighbourhood_fct=vshd_neighbourhood,
             **kwargs
):
    if len(hp_sequence) != len(hp_coordinates):
        raise ValueError("hp_sequence and hp_coordinates "
                         "needs to have the same length.")

    coord_i = np.copy(hp_coordinates)
    n = len(hp_sequence)
    for i in range(steps):
        selected_res_idx = random.randint(0, n-1)
        changed, coord_iplus1 = mc_move(coord_i, selected_res_idx, neighbourhood_fct)
        c_change_prob = metropolis_criterion(hp_sequence, coord_i, coord_iplus1, T, **kwargs)
        if changed:
            if c_change_prob == 1:
                coord_i = coord_iplus1
            elif random.random() <= c_change_prob:
                coord_i = coord_iplus1

    return coord_i


def REMCSimulation(conformations, optimal_energy, max_iter,
                   steps, neighbourhood_fct=vshd_neighbourhood
):  
    if True:
        sorted_conformation = sorted(conformations, key=lambda x: x.T)
        if not (conformations == sorted_conformation):
            raise ValueError(f"Temperature values are not sorted: {[i.T for i in conformations]}")

    hp_sequence = conformations[0].hp_sequence
    n_iter = 0
    best_energy = 0
    offset = 0
    while((best_energy > optimal_energy) & (n_iter < max_iter)):
        n_iter += 1
        for replica in conformations:
            replica.search(steps=steps, neighbourhood_fct=vshd_neighbourhood)
            if replica.energy < best_energy:
                best_energy = replica.energy
        i = offset
        while(i+1 < len(conformations)):
            j = i + 1
            re_probability = re_criterion(hp_sequence=hp_sequence,
                                          hp_coordinates_i=conformation[i], T_i=conformation[i].T,
                                          hp_coordinates_j=conformation[j], T_j=conformation[j].T,
                                          **kwargs)
            if re_probability == 1:
                conformation[i].swapTemperature(conformation[j])
            elif random.random() <= re_probability:
                conformation[i].swapTemperature(conformation[j])
            i = i + 2
        offset = 1 - offset

    return conformations


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
