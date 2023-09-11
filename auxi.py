H_RES = ("V", "I", "F", "L", "M", "C", "W")
OTHER_RES = ("D", "E", "K", "R", "H", "Y", "S", "T", "N", "Q",
            "G", "A", "P",
            "B", "Z", "X", "J", "O", "U")  # specials

HP_RES = {
    **dict.fromkeys(H_RES, "H"),
    **dict.fromkeys(OTHER_RES, "P")
    }

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


def is_H(converted_residue):
    return converted_residue == "H"


def random_index(positions):
    """Return a random index from the available indices."""
    selected_index = random.randint(0, len(positions) - 1)
    return selected_index
