import random
import itertools
from protocellule.specie import specie
from protocellule.reazione import reazione


class ChemGenerator:

    def __init__(self, conf_file="generator_conf.txt"):
        self.conf_file = conf_file
        self.default_buffered_specie = dict()

        self.specie_pool = list()
        self.catalyst_pool = list()
        self.reactions = list()

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
                self.default_buffered_specie[line.split()[0]] = float(line.split()[1])

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

        symbol_list = list(map(ChemGenerator.number_to_letter, combs))

        self.specie_pool = self.letter_to_specie(symbol_list)
        self.fill_catalyst()    # side-effect on self.catalyst_pool

        # for el in self.specie_pool:
        #     print(el)

    def generate_reactions(self):

        # generate membrane passage reaction
        for el in self.specie_pool:
            if el.inter is not None and el.inter[1]:
                (reagenti, prodotti, costante) = self.make_reaction(210, el)
                reaction = [reagenti, prodotti, costante]
                self.reactions.append(reazione(210, reaction))

        for i in range(self.n_cleavage):
            (reagenti, prodotti, costante) = self.make_reaction(23)
            reaction = list(reagenti + prodotti)
            reaction.append(costante)
            self.reactions.append(reazione(23, reaction))

        for i in range(self.n_condensation):
            (reagenti, prodotti, costante) = self.make_reaction(32)
            reaction = list(reagenti + prodotti)
            reaction.append(costante)
            self.reactions.append(reazione(32, reaction))

        # print(len(self.reactions))

    def generate_file(self):
        pass

    def generate_chem(self):
        self.generate_alphabet()
        self.generate_reactions()
        self.generate_file()

    def letter_to_specie(self, symbol_list):
        """Convert all symbls to the corresponding specie using the 'specie' Class
        The function also adds the buffered specie based on the 'default_buffered_specie'
        """
        out = list()
        for s in symbol_list:
            if s not in self.default_buffered_specie.keys():
                out.append(specie(s, self.init_specie_concentration))
            else:
                # add specie passing the membrane and the corresponding buffered external specie
                out.append(specie(s, self.init_specie_concentration, (0, 1)))
                out.append(specie(s+"ext", self.default_buffered_specie[s], (1, 0)))
        return out

    # works by side-effect on self.catalyst_pool
    def fill_catalyst(self):
        """Fill the catalyst_pool by looking at specie_pool and catalyst_min/max_length"""
        for el in self.specie_pool:
            if self.catalyst_min_length <= len(el.nome) <= self.catalyst_max_length:
                self.catalyst_pool.append(el)

    def make_reaction(self, reaction_type, el=None):
        """Generate a random reaction
        Returns a tuple of 3 element:
            'reagent': a tuple containing the names of the reagents
            'products': a tuple containing the names of the products
            'constant': the reaction constant
        """
        if reaction_type == 23 or reaction_type == 32:
            s = random.choice([el for el in self.specie_pool if len(el.nome) >= 2]).nome
            catalyst = random.choice(self.catalyst_pool).nome

            cut_point = random.randint(1, len(s)-1)

            first_part = s[0:cut_point]
            second_part = s[cut_point:len(s)]

            # cleavage
            if reaction_type == 23:
                return (s, catalyst), (first_part, second_part, catalyst), self.coeff_cleavage

            # condensation
            if reaction_type == 32:
                return (first_part, second_part, catalyst), (s, catalyst), self.coeff_condensation

        # membrane passage
        if reaction_type == 210 and el is not None:
            reagent = self.specie_pool[self.specie_pool.index(el.nome+"ext")]
            return reagent.nome, el.nome, self.coeff_membrane

    @staticmethod
    def number_to_letter(el):
        """Convert a number to the corresponding letter (0 -> A, 1 -> B, 2 -> C, ...)"""
        out = ""
        for i in range(len(el)):
            out += (chr(ord('A') + el[i]))
        return out


if __name__ == "__main__":
    # file = input("Inserisci nome file di configurazione (file.txt): ")
    # gen = ChemGenerator(file)
    gen = ChemGenerator()
