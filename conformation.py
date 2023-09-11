import numpy as np
from auxiliaries import *
from moves import *

class Conformation:
    def __init__(self, sequence, T=1, **kwargs):
        self.kwargs = kwargs
        name = kwargs.get("name", "unknown")
        hp_coordinates = kwargs.get("hp_coordinates", None)
        random_coord_initialization = kwargs.get("random", False)
        dtype = kwargs.get("dtype", np.int16)

        if hp_coordinates is not None:
            if len(sequence) != len(hp_coordinates):
                raise ValueError("Provided hp_coordinates arguments "
                                 "does not match sequence length")
        
        self.name = name
        self._sequence = sequence
        self._T_ini = T  # initial temperature
        self._T = T  # temperature
        self._hp_sequence = sequence_to_HP(self._sequence)
        self._hp_coordinates = hp_coordinates if hp_coordinates is not None else \
            initialize_coordinates(self._hp_sequence,
                                   random=random_coord_initialization,
                                   dtype=dtype)
        self._energy = conformation_energy(self._hp_sequence, self._hp_coordinates)

    def __str__(self):
        to_return = f"seq={self.sequence}, T={self.T}, energy={self.energy}\n" \
                    f"hp_seq={self.hp_sequence}\n" \
                    f"coord={self.hp_coordinates.tolist()}\n"
        return to_return

    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self, seq):
        # The sequence should not be modified this way
        self._sequence = seq
        self._hp_sequence = sequence_to_HP(self._sequence)
        self._energy = conformation_energy(self._sequence, self._hp_sequence)

    @property
    def hp_sequence(self):
        return self._hp_sequence

    @hp_sequence.setter
    def hp_sequence(self, value):
        # The hp_sequence should not be modified this way
        pass

    @property
    def hp_coordinates(self):
        return self._hp_coordinates

    @hp_coordinates.setter
    def hp_coordinates(self, value):
        self._hp_coordinates = value
        print(f"energy_i-1 = {self._energy}")
        self._energy = conformation_energy(self._hp_sequence, self._hp_coordinates, **self.kwargs)
        print(f"energy_i+1 = {self._energy}")

    @property
    def T(self):
        return self._T

    @T.setter
    def T(self, temperature):
        self._T = temperature

    @property
    def energy(self):
        return self._energy

    @energy.setter
    def energy(self, value):
        """Cannot change energy without changing conformation."""
        pass

    # Main functions
    def initial_temperature(self):
        return self._T_ini

    def search(self, steps, neighbourhood_fct):
        self.hp_coordinates = MCsearch(hp_sequence=hp_sequence,
                                       hp_coordinates=replica.hp_coordinates,
                                       T=replica.T,
                                       steps=steps,
                                       neighbourhood_fct=neighbourhood_fct,
                                       **self.kwargs)
        return self.hp_coordinates

    def swapTemperature(self, conformation2):
        T1, T2 = self.initial_temperature(), conformation2.initial_temperature()
        if T1 != T2:
            self.T = T2
            conformation2.T = T1

    def replica_exchange(self, conformation2):
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
