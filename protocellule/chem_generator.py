import random
import itertools
from protocellule.specie import Specie
from protocellule.reazione import Reazione


class ChemGenerator:

    def __init__(self, conf_file="generator_conf.txt"):
        self.conf_file = conf_file
        self.default_buffered_specie = dict()

        self.specie_pool = list()
        self.catalyst_pool = list()
        self.reactions = list()

        self.symm_specie_pool = list()
        self.symm_reactions = list()

        with open(self.conf_file) as f:
            riga = f.readline()
            self.file_output = riga.split()[0]+".txt"

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
            self.bimolecular = bool(int(riga.split()[0]))

            riga = f.readline()
            self.symmetric = bool(int(riga.split()[0]))

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

    def generate_chem(self):
        self.generate_alphabet()
        self.generate_reactions()

        self.delete_species(self.specie_pool, self.reactions)
        self.delete_species(self.symm_specie_pool, self.symm_reactions)

        self.generate_file(self.file_output, self.specie_pool, self.reactions)
        self.generate_file("SYMM_"+self.file_output, self.symm_specie_pool, self.symm_reactions)

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

        self.specie_pool.sort()

        # for el in self.specie_pool:
        #     print(el)

    def generate_reactions(self):
        """Generate random reactions
        The type of reactions are:
            210: membrane passage
            23: cleavage
            32: condensation
        """
        # generate membrane passage reaction
        for s in self.specie_pool:
            if s.inter is not None and s.inter["bufferizzata"]:
                (reagenti, prodotti, costante) = self.make_reaction(210, s)
                reaction = [reagenti, prodotti, costante]
                self.reactions.append(Reazione(210, reaction))

        for i in range(self.n_cleavage):
            (reagenti, prodotti, costante) = self.make_reaction(23)
            reaction = list(reagenti + prodotti)
            reaction.append(costante)
            self.reactions.append(Reazione(23, reaction))

        if self.symmetric:
            self.symm_specie_pool = self.specie_pool.copy()
            self.symm_reactions = self.reactions.copy()

        for i in range(self.n_condensation):
            (reagenti, prodotti, costante) = self.make_reaction(32)
            reaction = list(reagenti + prodotti)
            reaction.append(costante)
            r = Reazione(32, reaction)

            if self.symmetric:
                (l_complex, l_reaction) = self.tri_to_duo(r, True)

                for c in l_complex:
                    self.symm_specie_pool.append(Specie(c, 0, (0, 0)))

                for l in l_reaction:
                    self.symm_reactions.append(l)

            if self.bimolecular:
                (complex, l_reaction) = self.tri_to_duo(r)

                # add complex specie to pool
                self.specie_pool.append(Specie(complex, 0, (0, 0)))

                for l in l_reaction:
                    self.reactions.append(l)
            else:
                self.reactions.append(Reazione(32, reaction))

    # Monkey patch on specie object
    def delete_species(self, specie_pool, reactions):
        """Delete species which not appears in any reaction
        Monkey patch the attribute 'useless' on the 'specie' object
        """
        for s in specie_pool:
            useless = True
            for r in reactions:
                if s.nome in r.prodotti or s.nome in r.reagenti:
                    useless = False
                    break
            s.useless = useless

    def generate_file(self, filename, specie_pool, reactions):
        with open(filename, "w") as f:
            f.write(f"num_specie\t{sum([1 for s in specie_pool if not s.useless])}\n")
            f.write(f"num_reac\t{len(reactions)}\n")
            f.write(f"memb_thick\t{self.membrane_thickness}\n")
            f.write(f"cont_rad\t{self.radius}\n")
            f.write(f"density \t{self.density}\n")
            f.write("\n")

            for s in specie_pool:
                if not s.useless:
                    if s.inter is not None:
                        inter = str([v for k, v in s.inter.items()]).replace(",", "\t").replace("[", "").replace("]", "")
                    else:
                        inter = "0\t0"
                    f.write(f"{s.nome}\t{s.qnt}\t{inter}\n")

            f.write("\n")

            for r in reactions:
                f.write(str(r.tipo) + "\t" +
                        str(r.reagenti).replace(",", " +").replace("[", "").replace("]", "").replace("'", "") + " > " +
                        str(r.prodotti).replace(",", " +").replace("[", "").replace("]", "").replace("'", "") + " ; " +
                        str(r.costante) + "\n")

    def letter_to_specie(self, symbol_list):
        """Convert all symbols to the corresponding specie using the 'specie' Class
        The function also adds the buffered specie based on 'default_buffered_specie'
        """
        out = list()
        for s in symbol_list:
            if s not in self.default_buffered_specie.keys():
                out.append(Specie(s, self.init_specie_concentration, (0, 0.01)))
            else:
                # add specie passing the membrane and the corresponding buffered external specie
                out.append(Specie(s, self.init_specie_concentration, (0, 0.01)))
                out.append(Specie(s+"ext", self.default_buffered_specie[s], (1, 0)))
        return out

    # works by side-effect on self.catalyst_pool
    def fill_catalyst(self):
        """Fill the catalyst_pool by looking at specie_pool and catalyst_min/max_length"""
        for el in self.specie_pool:
            if (self.catalyst_min_length <= len(el.nome) <= self.catalyst_max_length) and "ext" not in el.nome and "*" not in el.nome:
                self.catalyst_pool.append(el)

    def make_reaction(self, reaction_type, buff_s=None):
        """Generate a random reaction
        Parameters:
            'reaction_type': number that specify the type of reaction
            'buff_s': specie object that is needed only to generate the membrane passage reaction
                      this is always a 'buffered' specie (the ones ending with 'ext')
        Returns a tuple of 3 element:
            'reagent': a tuple containing the names of the reagents
            'products': a tuple containing the names of the products
            'constant': the reaction constant
        """
        if reaction_type == 23 or reaction_type == 32:
            s = random.choice([el for el in self.specie_pool if len(el.nome) >= 2 and "ext" not in el.nome and "*" not in el.nome]).nome
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
        if reaction_type == 210 and buff_s is not None:
            prod_name = ''.join(buff_s.nome.rsplit('ext', 1))
            prod = self.specie_pool[self.specie_pool.index(prod_name)]
            return buff_s.nome, prod.nome, self.coeff_membrane

    def tri_to_duo(self, tri_reaction, symm=False):
        """Return the duo corresponding reaction of a tri reaction and the complex
        This method make sense only for condensation reaction and it returns a list with three reaction object
        Parameters:
            tri_reaction: the tri-reaction list
            symm: True if it's needed to generate the symmetric reactions
        """
        if tri_reaction.tipo != 32:
            return tri_reaction

        catalyst = tri_reaction.reagenti[-1]
        first_substr = tri_reaction.reagenti[0]
        second_substr = tri_reaction.reagenti[1]

        complex = '*' + first_substr + catalyst
        while complex in self.specie_pool:
            complex = '*' + complex

        first_r = Reazione(21, (first_substr, catalyst, complex, self.coeff_condensation/10))
        second_r = Reazione(12, (complex, first_substr, catalyst, 14.86))
        third_r = Reazione(22, (complex, second_substr, tri_reaction.prodotti[0], catalyst, self.coeff_condensation/20))

        if symm:
            f_complex = '*' + first_substr + catalyst
            while f_complex in self.symm_specie_pool:
                f_complex = '*' + f_complex

            s_complex = '*' + second_substr + catalyst
            while s_complex in self.symm_specie_pool:
                s_complex = '*' + s_complex

            first_r = Reazione(21, (first_substr, catalyst, f_complex, self.coeff_condensation / 10))
            second_r = Reazione(12, (f_complex, first_substr, catalyst, 14.86))
            third_r = Reazione(22, (f_complex, second_substr, tri_reaction.prodotti[0], catalyst, self.coeff_condensation / 20))

            symm_first_r = Reazione(21, (second_substr, catalyst, s_complex, self.coeff_condensation / 10))
            symm_second_r = Reazione(12, (s_complex, second_substr, catalyst, 14.86))
            symm_third_r = Reazione(22, (s_complex, first_substr, tri_reaction.prodotti[0], catalyst, self.coeff_condensation / 20))

            return (f_complex, s_complex), (first_r, second_r, third_r, symm_first_r, symm_second_r, symm_third_r)
        else:
            return complex, (first_r, second_r, third_r)

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
