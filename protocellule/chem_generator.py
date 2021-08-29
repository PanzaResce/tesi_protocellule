import random


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
        pass

    def generate_reactions(self):
        pass

    def generate_file(self):
        pass

    def generate_chem(self):
        self.generate_alphabet()
        self.generate_reactions()
        self.generate_file()


if __name__ == "__main__":
    # file = input("Inserisci nome file di configurazione (file.txt): ")
    # gen = ChemGenerator(file)
    gen = ChemGenerator()
