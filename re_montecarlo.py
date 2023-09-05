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
    """Read a fasta sequence and returns its sequence."""
    sequence = ""
    with open(filename, "r") as fasta_in:
        header = fasta_in.readline().strip()
        sequence = sequence + fasta_in.read().strip()
        sequence = header + sequence if header[0] != ">" else sequence
        print(sequence)
    return sequence


def seq_to_hp(sequence):
    """Convert protein sequence to its corresponding HP sequence."""
    pass


if __name__ == "__main__":
    filename = "./data/A0A0C5B5G6.fasta"
    read_fasta(filename)