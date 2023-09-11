
class Conformation:
    def __init__(sequence, T=1, **kwargs):
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
        self.__T = T  # initial temperature
        self._T = T  # temperature
        self._hp_sequence = sequence_to_HP(self._sequence)
        self._hp_coordinates = hp_coordinates if hp_coordinates is not None else \
            initialize_coordinates(self._hp_sequence,
                                   random=random_coord_initialization,
                                   dtype=dtype)
        self._energy = conformation_energy(self._hp_sequence, self._hp_coordinates)

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
        self._energy = conformation_energy(self._hp_sequence, self._hp_coordinates, **self.kwargs)

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

    def search(self, steps, neighbourhood_fct):
        self.hp_coordinates = MCsearch(hp_sequence=hp_sequence,
                                       hp_coordinates=replica.hp_coordinates,
                                       T=replica.T,
                                       steps=steps,
                                       neighbourhood_fct=neighbourhood_fct,
                                       **self.kwargs)
    return self.hp_coordinates