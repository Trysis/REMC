
class Conformation:
    def __init__(sequence, T=1, **kwargs):
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
        self._hp_coordinates = \
            initialize_coordinates(self._hp_sequence,
                                   random=random_coord_initialization,
                                   dtype=dtype)


        
    @property
    def sequence(self):
        return self._sequence

    @sequence.setter
    def sequence(self):
        # The sequence should not be modified this way
        pass

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
        pass
    
    @property
    def T(self, initial_temp=False):
        if initial_temp:
            return self.__T
        return self._T

    @T.setter
    def T(self, temperature):
        self._T = temperature
