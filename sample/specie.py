class specie:

    def __init__(self, nome, q, inter=None):
        ''' When inter = None it means that the species
        is involved in a reaction '''

        self.nome = nome
        self.qnt = float(q)

        if inter is not None:
            self.inter = {
                "bufferizzata": False,
                "boundary": False
            }

            if len(inter) != len(self.inter.keys()):
                raise LengthError(f"The interaction must be {len(self.inter.keys())}")

            if isinstance(inter, list):
                for index, key in enumerate(self.inter.keys()):
                    self.inter[key] = bool(int(inter[index]))
            else:
                self.inter = inter

        else:
            self.inter = None

    def __eq__(self, other):
        if isinstance(other, str):
            return self.nome == other
        else:
            return self.nome == other.nome and self.inter == other.inter

    def __str__(self):
        return f"\nNome specie: {self.nome}\nQuantità: {self.qnt}\nInterazioni: {self.inter}\n"

    def __repr__(self):
        return self.__str__()
