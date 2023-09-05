""""""

H_RES = ["V", "I", "F", "L", "M", "C", "W"]
OTHER_RES = ["D", "E", "K", "R", "H", "Y", "S", "T", "N", "Q"] + \
            ["G", "A", "P"] + \
            ["B", "Z", "X", "J", "O", "U"]  # specials

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


def seq_to_hp(sequence):
    """Convert protein sequence to its corresponding HP sequence."""
    sequence = sequence.upper()
    hp_sequence = ""
    for res in sequence:
        hp_sequence += HP_RES[res]
    return hp_sequence


def available_adjacentPos(hp_coordinates, i):
    """Returns free position adjacent to residue i."""
    

def init_coordinates(hp_sequence, random=False):
    positions = [[0, 0]]  # position of first residue
    if not random:
        for i in range(1, len(hp_sequence)):
            positions.append([i, 0])
    else:
        pass
    return np.array(positions)

if __name__ == "__main__":
    filename = "./data/A0A0C5B5G6.fasta"
    sequence = read_fasta(filename)
    hp_sequence = seq_to_hp(sequence)
    hp_coordinates = init_coordinates(hp_sequence)

    print(hp_coordinates)