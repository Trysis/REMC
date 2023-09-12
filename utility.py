"""This class contains auxiliaries function useful for multiple cases."""

# Residue 1 letter code associated with 3 letter code
RES_TO_3CODE = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
            'PYL': 'O', 'SEC': 'U'}  # Without B, Z, or unknown X one code letters

RES_TO_1CODE = {value: key for key, value in RES_TO_3CODE.items()}

H_RES = ("V", "I", "F", "L", "M", "C", "W")
OTHER_RES = ("D", "E", "K", "R", "H", "Y", "S", "T", "N", "Q",
            "G", "A", "P",
            "B", "Z", "X", "J", "O", "U")  # specials

HP_RES = {
    **dict.fromkeys(H_RES, "H"),
    **dict.fromkeys(OTHER_RES, "P")
    }


def read_fasta(filename):
    """Read the first fasta sequence and returns its sequence.

    filename: str
        Path to a fasta file (or valid sequence)

    Returns: str
        Return the associated sequence.

    """
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


def is_H(converted_residue):
    """Check if the argument is hydrophobic."""
    return converted_residue == "H"


def sequence_to_HP(sequence):
    """Convert protein sequence to its corresponding HP sequence."""
    sequence = sequence.upper()
    hp_sequence = ""
    for res in sequence:
        hp_sequence += HP_RES[res]
    return hp_sequence


def atom_section(atom_id, atom_name, res_1_name, res_id, x, y, z,
                 chain="A", occupancy=0, bfactor=0
):
    """Return a correct atom section format.

    atom_id: int (or str -> representing a numerical value)
        Atom id

    atom_name: str
        The given atom name

    res_1_name: str
        One letter residue name

    res_id: int (or str -> representing a numerical value)
        Residue id

    x, y, z: float (or int)
        x, y and z coordinate.

    chain: str
        Associated chain

    occupancy: int
        occupancy

    bfactor: int
        bfactor

    Returns: str
        a record ATOM line.

    """
    atom_type = atom_name[0]
    res_3_name = RES_TO_1CODE[res_1_name]
    atm_record = "" # Contiendra une section ATOM complète
    atm_record += f"{'ATOM':6s}{int(atom_id):5d} {atom_name:^4s}{'':1s}{res_3_name:3s} {chain:1s}"
    atm_record += f"{int(res_id):4d}{'':1s}   {x:8.3f}{y:8.3f}{z:8.3f}{float(occupancy):6.2f}"
    atm_record += f"{float(bfactor):6.2f}          {atom_type:>2s}{'':2s}"

    return atm_record


if __name__ == "__main__":
    atm_section = atom_section(atom_id=1, atom_name="CA",
                               res_1_name="G", res_id=1,
                               x=0.0, y=0.0, z=0.0,
                               bfactor=1.0)
    print(atm_section)