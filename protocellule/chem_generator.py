import random
import itertools


class ChemGenerator:

    def __init__(self, conf_file="generator_conf.txt"):
        self.conf_file = conf_file
        self.default_specie_concentration = dict()

        self.specie_pool = list()
        self.catalyst_pool = list()

        with open(self.conf_file) as f:
            riga = f.readline()
            self.seed = int(riga.split()[0])
            random.seed(self.seed)

            riga = f.readline()
            self.specie_max_length = int(riga.split()[0])

            riga = f.readline()
            self.alphabet_card = int(riga.split()[0])

            riga = f.readline()
            self.catalyst_min_length = int(riga.split()[0])

            riga = f.readline()
            self.catalyst_max_length = int(riga.split()[0])

            riga = f.readline()
            self.n_condensation = int(riga.split()[0])

            riga = f.readline()
            self.n_cleavage = int(riga.split()[0])

            riga = f.readline()
            self.init_specie_concentration = float(riga.split()[0])

            riga = f.readline()
            self.coeff_cleavage = int(riga.split()[0])

            riga = f.readline()
            self.coeff_condensation = int(riga.split()[0])

            riga = f.readline()
            self.coeff_membrane = float(riga.split()[0])

            riga = f.readline()
            self.membrane_thickness = float(riga.split()[0])

            riga = f.readline()
            self.radius = float(riga.split()[0])

            riga = f.readline()
            self.density = float(riga.split()[0])

            f.readline()    # empty row
            f.readline()    # title row

            for line in f:
                self.default_specie_concentration[line.split()[0]] = float(line.split()[1])

        # chem generation
        self.generate_chem()

    def generate_alphabet(self):
        """Generate all possible combinations of length L of X symbol
        L = self.specie_max_length
        X = self.alphabet_card
        The symbols are later converted from numbers to letters
        """
        combs = set()
        symbols = list(range(self.alphabet_card))

        for sub_length in range(1, self.specie_max_length+1):
            output = itertools.combinations_with_replacement(symbols, sub_length)

            for vett in output:
                combs.add(vett)
                all_combs = itertools.permutations(vett)
                for el in all_combs:
                    combs.add(el)

        self.specie_pool = list(map(ChemGenerator.number_to_letter, combs))

        print(len(self.specie_pool))

        for el in self.specie_pool:
            print(el)

    def generate_reactions(self):
        pass

    def generate_file(self):
        pass

    def generate_chem(self):
        self.generate_alphabet()
        self.generate_reactions()
        self.generate_file()

    @staticmethod
    def number_to_letter(el):
        """Convert a number to the corresponding letter (0 -> A, 1 -> B, 2 -> C, ...)"""
        out = list()
        for i in range(len(el)):
            out.append(chr(ord('A') + el[i]))
        return out


if __name__ == "__main__":
    # file = input("Inserisci nome file di configurazione (file.txt): ")
    # gen = ChemGenerator(file)
    gen = ChemGenerator()
