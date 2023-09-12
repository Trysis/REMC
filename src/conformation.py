"""Conformation representation."""

# Data gestion
import numpy as np

# Plot animation
from matplotlib import animation

# Local modules
from utility import *
from moves import *


class Conformation:
    """Class to represent a conformation."""

    def __init__(self, sequence, T=1, **kwargs):
        """Conformation representation.
        
        sequence: str
            The given protein sequence.

        T: int (or float)
            Given temperature for the conformation.
        
        **kwargs: dict, optional
            name: str
                Associated label to the conformation.

            hp_coordinates: numpy.ndarray -> of shape (n, 2)
                array-like containing x and y coordinates for each
                of the n residues.

            random: bool
                If random is False then the sequence is placed linearly,
                such that it begins from coordinates [0, 0], to [0, n-1]
                Else, it places residue randomly such that it begins from
                coordinates [0, 0] but end up at a random position [x, y].

            dtype: type
                numpy array type for the coordinates
            
            h_penalty: int (or float)
                Applied penalty when H residue are topological neighbour.

            o_penalty: int (or float)
                Applied score on the else statement of the penalty.

        Returns: None
            Instanciate the class with the given arguments.

        """
        # Argument retrieving
        name = kwargs.get("name", "unknown")
        hp_coordinates = kwargs.get("hp_coordinates", None)
        random_coord_initialization = kwargs.get("random", False)
        dtype = kwargs.get("dtype", np.int16)
        size = kwargs.get("size", (8, 8))

        if hp_coordinates is not None:
            if len(sequence) != len(hp_coordinates):
                raise ValueError("Provided hp_coordinates arguments "
                                 "does not match sequence length")

        if T <= 0:
            raise ValueError(f"Temperature={T} should be positive.")

        if not isinstance(sequence, str):
            raise ValueError(f"The given sequence {type(sequence) = } "
                             f"should be an str.")

        if not isinstance(name, str):
            raise ValueError(f"The given name {type(name) = } "
                             f"should be an str.")

        if not isinstance(random_coord_initialization, bool):
            raise ValueError(f"The given random type={type(name)} "
                             f"should be a bool.")

        self.kwargs = kwargs
        self.subplot = plt.subplots(1, 1, figsize=size)
        self.name = name
        self._sequence = sequence
        self.__T = T  # initial temperature
        self._T = T  # temperature
        self._hp_sequence = sequence_to_HP(self._sequence)
        self.__hp_coordinates = hp_coordinates if hp_coordinates is not None else \
            initialize_coordinates(hp_sequence=self._hp_sequence,
                                   random=random_coord_initialization,
                                   dtype=dtype)  # initial coordinates

        self._hp_coordinates = np.copy(self.__hp_coordinates)
        self._energy = conformation_energy(self._hp_sequence, self._hp_coordinates)
        self._transitions = []
        self.change_list = []
        self.T_list = [(0, self.T)]
        self.E_list = [self._energy]
        self._index = 1  # count on move
        self.iteration = 1  # count all iterations

    def __str__(self):
        """Redefine the output of the function when printed (with print function)."""
        to_return = f"______Conformation -{self.name}-______\n\n" \
                    f"sequence={self.sequence}\n" \
                    f"hp_sequence={self.hp_sequence}\n" \
                    f"temperature={self.T}\n" \
                    f"energy={self.energy}\n" \
                    f"coordinates={self.hp_coordinates.tolist()}\n" \
                    f"_____________\n"

        return to_return

    @property
    def index(self):
        return self._index

    @index.setter
    def index(self, value):
        pass

    @property
    def sequence(self):  # Getter
        return self._sequence

    @sequence.setter
    def sequence(self, seq):  # Setter
        # The sequence should not be modified this way
        self._sequence = seq
        self._hp_sequence = sequence_to_HP(self._sequence)
        self._energy = conformation_energy(self._sequence, self._hp_sequence)

    @property
    def hp_sequence(self):  # Getter
        return self._hp_sequence

    @hp_sequence.setter
    def hp_sequence(self, value):  # Setter
        # The hp_sequence should not be modified this way
        pass

    @property
    def hp_coordinates(self):  # Getter
        return self._hp_coordinates

    @hp_coordinates.setter
    def hp_coordinates(self, value):  # Setter
        self._hp_coordinates = value
        self._energy = conformation_energy(self._hp_sequence, self._hp_coordinates, **self.kwargs)

    @property
    def T(self):  # Getter
        return self._T

    @T.setter
    def T(self, temperature):  # Setter
        self._T = temperature

    @property
    def energy(self):  # Getter
        return self._energy

    @energy.setter
    def energy(self, value):  # Setter
        """Cannot change energy without changing conformation."""
        pass

    @property
    def transitions(self):
        return self._transitions

    @transitions.setter
    def transitions(self, value):
        if isinstance(value, np.ndarray):
            self._transition = value.tolist()
        if isinstance(value, list):
            self._transition = value

        raise ValueError("Incorrect type, transitions "
                         "should be an array.")

    # Main functions
    def initial_temperature(self):
        """Returns the initial temperature of the conformation."""
        return self.__T

    def search(self, steps, neighbourhood_fct, move_on_step):
        """Perform the Monte Carlo conformation search on the conformation.

        steps: int
            Number of steps for the Monte Carlo search algorithm.

        neighbourhood_fct: function
            The function needs to be an implemented VSHD, PULL_MOVE or other
            neighbourhood function. It should returns a dictionnary containing
            a set of available move for the chosen function. The dictionnary
            is structured such as dict[key] = list, list.
            The first list contains residue indices, and the second the coordinates.

        Returns: numpy.ndarray
            Returns the new coordinates after {steps} iteration of
            tried fold on the conformation {i}.

        """
        for i in range(steps):
            changed, self.hp_coordinates = MCsearch(hp_sequence=self.hp_sequence,
                                                    hp_coordinates=self.hp_coordinates,
                                                    T=self.T,
                                                    steps=1,
                                                    neighbourhood_fct=neighbourhood_fct,
                                                    move_on_step=move_on_step,
                                                    **self.kwargs)

            self.iteration += 1

            if changed:
                self._index += 1
                self._transitions.append(self._get_artist())
                self.change_list.append(self.iteration-1)
                self.E_list.append(self.energy)

        return self.hp_coordinates

    def swapTemperature(self, conformation2):
        """Exchange the temperature between two conformation.
        
        conformation2: Conformation
            A conformation to exchange temperature with.

        Returns: None

        """
        T1, T2 = self.T, conformation2.T
        if T1 != T2:
            self.T = T2
            conformation2.T = T1

    def replica_exchange(self, conformation2):
        """Exchange temperature between two conformation
            with a certain probability.

        conformation2: Conformation
            A conformation to exchange temperature with.

        Returns: None

        """
        re_probability = re_criterion(hp_sequence=self.hp_sequence,
                                      hp_coordinates_i=self.hp_coordinates,
                                      T_i=self.T,
                                      hp_coordinates_j=conformation2.hp_coordinates,
                                      T_j=conformation2.T,
                                      **self.kwargs)

        # Swap temperature with a certain probability
        if re_probability == 1:
            self.swapTemperature(conformation2)
            self.T_list.append((self.iteration, self.T))
        elif random.random() <= re_probability:
            self.swapTemperature(conformation2)
            self.T_list.append((self.iteration, self.T))

    def plot(self, **kwargs):
        """"""
        kwargs["T"] = kwargs.get("T", self.T)
        kwargs["legend_title"] = kwargs.get("legend_title", "Conformation")
        _, ax = plot_conformation(self.hp_coordinates, self.hp_sequence, show=True, **kwargs)

    def _get_artist(self):
        """"""
        all_artist = plot_conformation(self.hp_coordinates, self.hp_sequence, subplot=self.subplot,
                                       show=False, returns_artist=True, 
                                       T=self.T, legend_title="Conformation",
                                       index=self._index, energy=self.energy,
                                       **self.kwargs)

        return all_artist

    # TODO: 2 calls will fail
    def animate(self, fps=6, save=False, show=True, **kwargs):
        if len(self._transitions) > 0:
            plt_animation = animation.ArtistAnimation(self.subplot[0], self._transitions, repeat=False)
            if save:
                plt_animation.save(f"../img/conf_{self.name}.gif", fps=fps, dpi=300)

            if show:
                plt.show()

# TODO: Generate Pymol coordinate (animation or single frame)

if __name__ == "__main__":
    conf = Conformation("APKGGAYKVVVVVVVVVVVVAP", name="test", T=1, random=False)
    conf.search(steps=200, neighbourhood_fct=vshd_neighbourhood, move_on_step=True)
    conf.animate(save=True)
