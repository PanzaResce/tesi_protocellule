class Specie:

    def __init__(self, nome, q, inter=None):
        """When inter = None it means that the species
        is involved in a reaction """

        self.nome = nome
        self.qnt = float(q)

        if inter is not None:
            self.inter = {
                "bufferizzata": 0.0,
                "boundary": 0.0
            }

            if len(inter) != len(self.inter.keys()):
                raise AttributeError(f"Wrong format, the interactions must be {len(self.inter.keys())}")

            if isinstance(inter, list) or isinstance(inter, tuple):
                for index, key in enumerate(self.inter.keys()):
                    self.inter[key] = float(inter[index])
            else:
                self.inter = inter

        else:
            self.inter = None

    def __eq__(self, other):
        if isinstance(other, str):
            return self.nome == other
        else:
            return self.nome == other.nome and self.inter == other.inter

    def __lt__(self, other):
        self_buffered = (self.inter is not None and self.inter["bufferizzata"] == 1)
        other_buffered = (other.inter is not None and other.inter["bufferizzata"] == 1)
        both_buffered = (self_buffered and other_buffered)

        both_equal_length = len(self.nome) == len(other.nome)

        if both_buffered:
            if both_equal_length:
                return self.nome < other.nome               # lexicographic
            else:
                return len(self.nome) < len(other.nome)     # shorter first
        elif self_buffered:                                 # buffered first
            return True
        elif other_buffered:
            return False
        else:                                               # if are both not buffered
            if both_equal_length:
                return self.nome < other.nome               # lexicographic
            else:
                return len(self.nome) < len(other.nome)     # shorter first

    def __str__(self):
        return f"\nNome specie: {self.nome}\nQuantità: {self.qnt}\nInterazioni: {self.inter}\n"

    def __repr__(self):
        return self.__str__()
