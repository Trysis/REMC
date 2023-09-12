"""This class contains function to associated with
    coordinates and moves applied on conformation.

"""

# Python default modules
import random

# Data gestion
import numpy as np
import matplotlib.pyplot as plt

# Local modules
from utility import *

# Coordinates type
COORD_TYPE = np.int16

# 2D
XY_UP = np.array([0, 1], dtype=COORD_TYPE)
XY_RIGHT = np.array([1, 0], dtype=COORD_TYPE)
XY_LEFT = np.array([-1, 0], dtype=COORD_TYPE)
XY_DOWN = np.array([0, -1], dtype=COORD_TYPE)


def random_index(positions):
    """Return a random index from the available indices.
    
    positions: list or numpy.ndarray of shape (n, 2)
        array-like containing x and y coordinates for
        n residues.

    Returns: int
        Random indice in [0, n-1]
    """
    selected_index = random.randint(0, len(positions) - 1)
    return selected_index


def initialize_coordinates(hp_sequence, random=False, dtype=np.int16):
    """Returns a set of coordinates for the associated sequence.

    hp_sequence: str
        Sequence in HP_model format containing H for hydrophobic residue
        or P for polar residue. Such that set(hp_sequence) == {H, P}

    random: bool
        If random is False then the sequence is placed linearly,
        such that it begins from coordinates [0, 0], to [0, n-1]
        Else, it places residue randomly such that it begins from
        coordinates [0, 0] but end up at a random position [x, y].

    dtype: type
        numpy array type for the coordinates

    Returns: numpy.ndarray -> of shape (n, 2)
        Returns the given coordinates initialized position associated
        to the {hp_sequence}.

    """
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
    """Returns available crank shaft move positions.
        {i} should be comprised in [1, n-1] and
        A sequence of length 3 or less cannot perform this move.
    
    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    Returns: list, list -> of shape: (x, 1), (x, 2)
        Returns two array, the first contains the residues that are
        to be moved and the second contains the coordinates where residues
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
    """Returns available positions to perform pull move for residue {i}.
        {i} should be comprised in [1, n-1] because {i-1} & {i+1} residue
        are needed.

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    Returns: list, list -> of shape: (x, 1), (x, 2)
        Returns two array, the first contains the residues that are
        to be moved and the second contains the coordinates where residues
        in the first array can move.

    """
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
    """Returns a set of available VSHD moves for residue {i}.
        Up to end_moves, corner_moves and crank-shaft moves.

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    Returns: dict -> dict[key] = list, list
        All available VSHD moves such as end_moves, corner_moves,
        and crank-shaft moves are returned for the residue {i}. The
        length of the dict can vary in [0, 3] depending of the moves
        that can be performed on residue {i}.

    """
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

    return available_movements


def pull_move_neighbourhood(hp_coordinates, i):
    """Returns an available pull move for residue {i}.

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    Returns: dict -> dict[key] = list, list
        An available pull move is returned for residue {i}. The
        length of the dict can vary in [0, 1] depending on the moves
        that can be performed on residue {i}.

    """
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

    return available_movements

def hybrid_neighbourhood(p=0.5):
    """Returns a random function between the chosen
        neighbourhood function.

    p: float
        Probability to use the Pull moves algorithm, else VSHD
        is performed.

    Returns: function
        It returns a function that can perform a move on residue i
        that will returns a dictionnary containing a set of available move
        for the chosen function. The dictionnary is structured such as
        dict[key] = list, list.
        The first list contains residue indices, and the second the coordinates.

    """
    if p == 1:
        return pull_move_neighbourhood
    if random.random() < p:
        return pull_move_neighbourhood

    return vshd_neighbourhood

def conformation_energy(hp_sequence, hp_coordinates, dtype=np.int16, **kwargs):
    """Get the energy from the corresponding conformation.

    hp_sequence: str
        Sequence in HP_model format containing H for hydrophobic residue
        or P for polar residue. Such that set(hp_sequence) == {H, P}

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    dtype: type
        numpy array type for the coordinates
    
    **kwargs:
        h_penalty: int (or float)
            Applied penalty when H residue are topological neighbour.
        o_penalty: int (or float)
            Applied score on the else statement of the penalty.

    Returns: int (or float)
        Returns the total energy of the given conformation.

    """
    if len(hp_sequence) != len(hp_coordinates):
        raise ValueError("arg1 and arg2 needs to have the same length.")
    if not isinstance(hp_coordinates, np.ndarray):
        hp_coordinates = np.array(hp_coordinates, dtype=dtype)
    # Penalty score
    H_penalty = kwargs.get("h_penalty", -1)
    otherwise_penalty = kwargs.get("o_penalty", 0)
    # Calculate energy
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
    """Metropolis criterion to change conformation from conformation {i} to {i+1}.

    hp_sequence: str
        Sequence in HP_model format containing H for hydrophobic residue
        or P for polar residue. Such that set(hp_sequence) == {H, P}

    hp_coordinates_i: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    hp_coordinates_iplus1: numpy.ndarray -> of shape (n, 2)
        Same as {hp_coordinates_i} but correspond to the new conformation after
        a conformational change has been performed on {hp_coordinates_i}.

    **kwargs:
        h_penalty: int (or float)
            Applied penalty when H residue are topological neighbour.
        o_penalty: int (or float)
            Applied score on the else statement of the penalty.

    Returns: float
        Returns the probability for the conformation {i} to change
        to conformation {i+1} (iplus1).

    """
    energy_i = conformation_energy(hp_sequence, hp_coordinates_i, **kwargs)
    energy_iplus1 = conformation_energy(hp_sequence, hp_coordinates_iplus1, **kwargs)
    deltaEnergy = energy_iplus1 - energy_i
    if deltaEnergy <= 0:
        return 1.0
    else:
        return np.exp((-deltaEnergy/T)).item()


def re_criterion(hp_sequence, hp_coordinates_i, T_i, hp_coordinates_j, T_j, **kwargs):
    """Replica exchange probability between two conformations {i} and {j}.

    hp_sequence: str
        Sequence in HP_model format containing H for hydrophobic residue
        or P for polar residue. Such that set(hp_sequence) == {H, P}

    hp_coordinates_i, hp_coordinates_j: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    T_i, T_j: int (or float)
        Given temperature to conformation {i} and conformation {j}.

    **kwargs:
        h_penalty: int (or float)
            Applied penalty when H residue are topological neighbour.
        o_penalty: int (or float)
            Applied score on the else statement of the penalty.

    Returns: float
        Returns the probability for the exchange of temperature of
        conformation {i} with conformation {j}.

    """
    energy_i = conformation_energy(hp_sequence, hp_coordinates_i, **kwargs)
    energy_j = conformation_energy(hp_sequence, hp_coordinates_j, **kwargs)
    criterion = ((1/T_j) - (1/T_i)) * (energy_i - energy_j)
    if criterion <= 0:
        return 1.0
    else:
        return np.exp(-criterion).item()


def mc_move(hp_coordinates, i, neighbourhood_fct=vshd_neighbourhood):
    """Try a movement on the conformation for residue {i}.

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    i: int
        index of the selected residue to get available end-moves
        positions on.

    neighbourhood_fct: function
        The function needs to be an implemented VSHD, PULL_MOVE or other
        neighbourhood function. It should returns a dictionnary containing
        a set of available move for the chosen function. The dictionnary
        is structured such as dict[key] = list, list.
        The first list contains residue indices, and the second the coordinates.
    
    Returns: bool, numpy.ndarray
        The function returns as its first argument True if a movement have been
        performed on residue {i} when available, or False otherwise.
        The second argument remains the same coordinates as {hp_coordinates} if
        no moves have been performed, else it contains new coordinates.

    """
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
    """Try multiple fold on the given conformation.

    hp_sequence: str
        Sequence in HP_model format containing H for hydrophobic residue
        or P for polar residue. Such that set(hp_sequence) == {H, P}

    hp_coordinates: numpy.ndarray -> of shape (n, 2)
        array-like containing x and y coordinates for each
        of the n residues.

    T: int (or float)
        Given temperature to the conformation {i}.

    steps: int
        Number of steps for the Monte Carlo search algorithm.

    neighbourhood_fct: function
        The function needs to be an implemented VSHD, PULL_MOVE or other
        neighbourhood function. It should returns a dictionnary containing
        a set of available move for the chosen function. The dictionnary
        is structured such as dict[key] = list, list.
        The first list contains residue indices, and the second the coordinates.

    **kwargs:
        h_penalty: int (or float)
            Applied penalty when H residue are topological neighbour.
        o_penalty: int (or float)
            Applied score on the else statement of the penalty.

    Returns: numpy.ndarray
        Returns the new coordinates after {steps} iteration of
        tried fold on the conformation {i}.

    """
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
    """Perform the Replica Exchange Monte Carlo Algorithm.

    conformations: list -> list(Conformation)
        Contains a list of conformations.

    optimal_energy: int (or float)
        The corresponding energy we want to stop to.

    max_iter: int
        Max iteration before stopping the algorithm, 
        used in case if optimal energy is never found.

    steps: int
        Steps used for the MCSearch algorithm, used to potentially
        perform {steps} folds on the given conformations.

    neighbourhood_fct: function
        The function needs to be an implemented VSHD, PULL_MOVE or other
        neighbourhood function. It should returns a dictionnary containing
        a set of available move for the chosen function. The dictionnary
        is structured such as dict[key] = list, list.
        The first list contains residue indices, and the second the coordinates.

    Returns: list -> list(Conformation)
        It returns the set of modified conformations after {steps}*{max_iter}
        iterations.

    """
    if True:
        sorted_conformation = sorted(conformations, key=lambda x: x.T)
        if not (conformations == sorted_conformation):
            raise ValueError(f"Temperature values are not sorted: {[i.T for i in conformations]}")

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
            conformations[i].replica_exchange(conformations[j])
            i = i + 2
        offset = 1 - offset

    return conformations


def plot_conformation(hp_coordinates, hp_sequence=None):
    """Plot the given 2D conformation."""
    plt.scatter(hp_coordinates[:, 0], hp_coordinates[:, 1])
    if hp_sequence is not None:
        if(len(hp_coordinates) != len(hp_sequence)):
            raise ValueError("arg1 and arg2 needs to have the same length.")
        for i, hp_letter in enumerate(hp_sequence):
            plt.text(hp_coordinates[i, 0], hp_coordinates[i, 1], hp_letter)
    plt.plot(hp_coordinates[:, 0], hp_coordinates[:, 1])


if __name__ == "__main__":
    pass
