"""This class contains auxiliaries function useful for multiple cases."""

# Data Gestion
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# 3 letter code residue to -> 1 letter code residue
RES_TO_1CODE = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
            'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N',
            'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W',
            'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M',
            'PYL': 'O', 'SEC': 'U'}  # Without B, Z, or unknown X one code letters

# 1 letter code to -> 3 letter code
RES_TO_3CODE = {value: key for key, value in RES_TO_1CODE.items()}

# HP Residus
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
    res_3_name = RES_TO_3CODE[res_1_name]
    atm_record = "" # Contiendra une section ATOM complète
    atm_record += f"{'ATOM':6s}{int(atom_id):5d} {atom_name:^4s}{'':1s}{res_3_name:3s} {chain:1s}"
    atm_record += f"{int(res_id):4d}{'':1s}   {x:8.3f}{y:8.3f}{z:8.3f}{float(occupancy):6.2f}"
    atm_record += f"{float(bfactor):6.2f}          {atom_type:>2s}{'':2s}"

    return atm_record


def legend_patch(label, color="none"):
    """Returns the corresponding label with a Patch object for matplotlib legend purposes.
    
    label: str
        Character chain specifying the label
    color: str
        Corresponding color attribute for label

    Returns: tuple {str, matplotlib.patches.Patch}
        Specified label with Patch

    """
    return mpatches.Patch(color=color, label=label), label


def plot_conformation(hp_coordinates, hp_sequence=None, show=False,
                      subplot=None, returns_artist=False, **kwargs
):
    """Plot the given 2D conformation."""

    # Argument retrieving
    size = kwargs.get("size", (10, 10))
    title = kwargs.get("title", "")
    xlabel = kwargs.get("xlabel", "")
    ylabel = kwargs.get("ylabel", "")
    top = kwargs.get("top", None)
    bottom = kwargs.get("bottom", None)
    temperature = kwargs.get("T", None)
    energy = kwargs.get("energy", None)
    legend_title = kwargs.get("legend_title", None)

    all_artist = []
    # Figure
    fig, ax = plt.subplots(1, 1, figsize=size) if subplot is None else subplot
    if returns_artist:
        all_artist.append(ax.scatter(hp_coordinates[:, 0], hp_coordinates[:, 1], color="grey"))
    else:
        ax.scatter(hp_coordinates[:, 0], hp_coordinates[:, 1])

    # Plot parameters
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_ylim(top=top, bottom=bottom)

    if hp_sequence is not None:
        if(len(hp_coordinates) != len(hp_sequence)):
            raise ValueError("arg1 and arg2 needs to have the same length.")
        for i, hp_letter in enumerate(hp_sequence):
            if returns_artist:
                all_artist.append(ax.text(hp_coordinates[i, 0], hp_coordinates[i, 1], hp_letter))
            else:
                ax.text(hp_coordinates[i, 0], hp_coordinates[i, 1], hp_letter)

    if returns_artist:
        im, = ax.plot(hp_coordinates[:, 0], hp_coordinates[:, 1], color="orange")
        all_artist.append(im)
    else:
        ax.plot(hp_coordinates[:, 0], hp_coordinates[:, 1])

    # Main Legend
    handles_main, labels_main = ax.get_legend_handles_labels()
    if (handles_main != []) & (labels_main != []):
        main_legend = ax.legend(handles_main, labels_main, loc="upper left")
        ax.add_artist(main_legend)

    # Loss Legend
    if True:
        handles, labels = [], []
        if temperature is not None:
            temp_patch, temp_label = legend_patch(f"T = {temperature:3.1f}")
            handles.append(temp_patch)
            labels.append(temp_label)
        
        if energy is not None:
            index_patch, index_label = legend_patch(f"E = {energy:<4.1f}")
            handles.append(index_patch)
            labels.append(index_label)

        if len(handles) > 0:
            set_legend = ax.legend(handles, labels, loc="best",
                                   title=legend_title,
                                   handlelength=0, handletextpad=0)

            if returns_artist:
                all_artist.append(set_legend)

            ax.add_artist(set_legend)

    ax.autoscale()

    if show:
        plt.show()

    if returns_artist:
        return all_artist

    return fig, ax


if __name__ == "__main__":
    atm_section = atom_section(atom_id=1, atom_name="CA",
                               res_1_name="G", res_id=1,
                               x=0.0, y=0.0, z=0.0,
                               bfactor=1.0)

    print(atm_section)
