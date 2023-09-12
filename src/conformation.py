"""Conformation representation."""

# Data gestion
import numpy as np

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
        self.kwargs = kwargs
        name = kwargs.get("name", "unknown")
        hp_coordinates = kwargs.get("hp_coordinates", None)
        random_coord_initialization = kwargs.get("random", False)
        dtype = kwargs.get("dtype", np.int16)

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
        print(f"energy_i-1 = {self._energy}")
        self._energy = conformation_energy(self._hp_sequence, self._hp_coordinates, **self.kwargs)
        print(f"energy_i+1 = {self._energy}")

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

    # Main functions
    def initial_temperature(self):
        """Returns the initial temperature of the conformation."""
        return self.__T

    def search(self, steps, neighbourhood_fct):
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
        self.hp_coordinates = MCsearch(hp_sequence=hp_sequence,
                                       hp_coordinates=replica.hp_coordinates,
                                       T=replica.T,
                                       steps=steps,
                                       neighbourhood_fct=neighbourhood_fct,
                                       **self.kwargs)

        return self.hp_coordinates

    def swapTemperature(self, conformation2):
        """Exchange the temperature between two conformation.
        
        conformation2: Conformation
            A conformation to exchange temperature with.

        Returns: None

        """
        T1, T2 = self.initial_temperature(), conformation2.initial_temperature()
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
                               hp_coordinates_i=self.hp_coordinates, T_i=self.T,
                               hp_coordinates_j=conformation2.hp_coordinates, T_j=conformation2.T,
                               **self.kwargs)
        if re_probability == 1:
            self.swapTemperature(conformation2)
        elif random.random() <= re_probability:
            self.swapTemperature(conformation2)

if __name__ == "__main__":
    conf1 = Conformation("APKGGAYK")
    print(conf1)
