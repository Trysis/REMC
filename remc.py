""""""
__authors__ = "Roude JEAN MARIE"
__contact__ = ("roude.etu@gmail.com", "roude.bioinfo@gmail.com")
__date__ = "01-09-2023"
__version__ = "1.0.0"

# Local module
from conformation import Conformation
from utility import *

if __name__ == "__main__":
    filename = "./data/A0A0C5B5G6.fasta"
    sequence = read_fasta(filename)
    hp_sequence = sequence_to_HP(sequence)
    hp_coordinates = initialize_coordinates(hp_sequence, random=True)
    print(f"{conformation_energy(hp_sequence, hp_coordinates) = }")
    plot_conformation(hp_coordinates, hp_sequence)
    #plt.show()
    hp_coo_2 = MCsearch(hp_sequence, hp_coordinates, T=1, steps=10000,
                        neighbourhood_fct=pull_move_neighbourhood)

    #plot_conformation(hp_coo_2, hp_sequence)
    #plt.show()
    conf1 = Conformation(sequence, T=1)
    conf2 = Conformatin(sequence, T=1)
    #REMCSimulation()
