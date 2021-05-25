import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from math import pi
from protocellule.specie import specie
from protocellule.reazione import reazione


class Proto:
    """ Per come è concepita la classe proto, tutte
        le informazioni che vengono memorizzate rappresentanto
        lo stato iniziale della cellula """

    def __init__(self, conf_file="conf.txt", chem_file="chem.txt"):

        self.specie = list()
        self.reazioni = list()
        self.history = dict()
        self.division_history = dict()
        # self.flow_history = list()

        self.conf_file = conf_file
        self.chem_file = chem_file

        with open(chem_file) as f:
            # Lettura proprietà cellula
            riga = f.readline()
            self.n_specie = int(riga.split()[1])

            riga = f.readline()
            self.n_reazioni = int(riga.split()[1])

            riga = f.readline()
            self.membrane_thickness = float(riga.split()[1])

            riga = f.readline()
            radius = float(riga.split()[1])

            # volume protocellula in litri
            self.volume = 4/3*pi*pow(radius, 3) * 1000

            riga = f.readline()
            self.density = float(riga.split()[1])

            # quantità iniziale lipide (mol) = densità * volume membrana
            self.contenitore = self.density * ((4/3)*pi*pow(radius + self.membrane_thickness, 3) - self.volume*0.001)

            f.readline()    #riga vuota

            # Lettura specie --> nome, quantità, proprietà
            for i in range(self.n_specie):
                riga = f.readline()
                s = specie(riga.split()[0], riga.split()[1], riga.split()[2:])

                self.specie.append(s)

            f.readline()   #riga vuota

            # Lettura reazioni
            for i in range(self.n_reazioni):
                riga = f.readline()

                tipo = riga.split()[0]
                vett_reazione = [el for el in riga.split()[1:] if el not in ("+", ">", ";", "/")]

                # this is for optimization purpose, the part (36pi^1/3 / membrane) is constant through all the program
                if tipo == "210":
                    vett_reazione[2] = float(vett_reazione[2]) * pow(36 * pi, 1 / 3) / self.membrane_thickness

                self.reazioni.append(reazione(tipo, vett_reazione))

    def simula(self):
        (x0, xn, n, n_div, div_coeff, print_div, print_hist) = Proto.read_conf(self.conf_file)
        self.div_coeff = div_coeff

        t_span = [x0, xn]

        y0 = [s.qnt for s in self.specie]
        y0.append(self.contenitore)

        for i in range(n_div):
            sol1 = solve_ivp(self.fn, t_span, y0, events=self.terminate)

            n_step = sol1.y.shape[1]
            out = [s[n_step - 1] for s in sol1.y]

            print(f"N. divisione: {i+1}")

            if print_div:
                self.fill_division_history(sol1.t[n_step-1], out)

            new_qnt = self.duplicate(out[-1])

            out[-1] = new_qnt
            y0 = out

            # va fatto un ulteriore controllo sulle quantità negative
            for idx in range(len(y0)):
                if y0[idx] < 0:
                    y0[idx] = 0

            if print_hist:
                self.fill_full_history(sol1.t, sol1.y)

        print("END")

    def fn(self, t, specie):
        """ Restituisce un vettore che indica le nuove
            quantità a seguito delle reazioni """

        # Calcolo volume attuale
        self.volume = self.calc_volume(specie[-1])

        # calcolo variazione quantità di lipide con vecchio volume
        dC = sum([s.inter["boundary"] * specie[self.specie.index(s)] for s in self.specie]) * self.volume

        delta = [0] * len(specie)

        if t != 0:
            # calcolo volume attuale con nuova quantità lipide
            new_volume = self.calc_volume(specie[-1] + dC)
            v_rapp = self.volume / new_volume

            # ricalcolo concentrazioni sostanze con rapporto (volume vecchio / volume nuovo)
            for idx, s in enumerate(specie):
                if idx != len(specie)-1 and bool(self.specie[idx].inter["bufferizzata"]) is not True:
                    delta[idx] -= specie[idx]*(1-v_rapp)

            self.volume = new_volume

        # applico le reazioni
        for i in range(self.n_reazioni):
            if self.reazioni[i].tipo == 12:
                i_reagente_0 = self.specie.index(self.reazioni[i].reagenti[0])
                i_prodotto_0 = self.specie.index(self.reazioni[i].prodotti[0])
                i_prodotto_1 = self.specie.index(self.reazioni[i].prodotti[1])
                flusso = self.reazioni[i].costante * specie[i_reagente_0]
                if specie[i_reagente_0] <= 0:
                    flusso = 0

                delta[i_reagente_0] -= flusso
                delta[i_prodotto_0] += flusso
                delta[i_prodotto_1] += flusso
            elif self.reazioni[i].tipo == 21:
                i_reagente_0 = self.specie.index(self.reazioni[i].reagenti[0])
                i_reagente_1 = self.specie.index(self.reazioni[i].reagenti[1])
                i_prodotto_0 = self.specie.index(self.reazioni[i].prodotti[0])
                flusso = self.reazioni[i].costante * specie[i_reagente_0] * specie[i_reagente_1]
                if specie[i_reagente_0] <= 0:
                    flusso = 0
                if specie[i_reagente_1] <= 0:
                    flusso = 0

                delta[i_reagente_0] -= flusso
                delta[i_reagente_1] -= flusso
                delta[i_prodotto_0] += flusso
            elif self.reazioni[i].tipo == 22:
                i_reagente_0 = self.specie.index(self.reazioni[i].reagenti[0])
                i_reagente_1 = self.specie.index(self.reazioni[i].reagenti[1])
                i_prodotto_0 = self.specie.index(self.reazioni[i].prodotti[0])
                i_prodotto_1 = self.specie.index(self.reazioni[i].prodotti[1])
                flusso = self.reazioni[i].costante * specie[i_reagente_0] * specie[i_reagente_1]
                if specie[i_reagente_0] <= 0:
                    flusso = 0
                if specie[i_reagente_1] <= 0:
                    flusso = 0

                delta[i_reagente_0] -= flusso
                delta[i_reagente_1] -= flusso
                delta[i_prodotto_0] += flusso
                delta[i_prodotto_1] += flusso
            elif self.reazioni[i].tipo == 10:
                i_reagente_0 = self.specie.index(self.reazioni[i].reagenti[0])
                flusso = self.reazioni[i].costante * specie[i_reagente_0]
                if specie[i_reagente_0] <= 0:
                    flusso = 0

                delta[i_reagente_0] -= flusso
            elif self.reazioni[i].tipo == 1:
                i_prodotto_0 = self.specie.index(self.reazioni[i].prodotti[0])
                flusso = self.reazioni[i].costante

                delta[i_prodotto_0] += flusso
            elif self.reazioni[i].tipo == 23:
                i_reagente_0 = self.specie.index(self.reazioni[i].reagenti[0])
                i_reagente_1 = self.specie.index(self.reazioni[i].reagenti[1])
                i_prodotto_0 = self.specie.index(self.reazioni[i].prodotti[0])
                i_prodotto_1 = self.specie.index(self.reazioni[i].prodotti[1])
                i_prodotto_2 = self.specie.index(self.reazioni[i].prodotti[2])
                flusso = self.reazioni[i].costante * specie[i_reagente_0] * specie[i_reagente_1]
                if specie[i_reagente_0] <= 0:
                    flusso = 0
                if specie[i_reagente_1] <= 0:
                    flusso = 0

                delta[i_reagente_1] -= flusso
                delta[i_reagente_1] -= flusso
                delta[i_prodotto_0] += flusso
                delta[i_prodotto_1] += flusso
                delta[i_prodotto_2] += flusso
            elif self.reazioni[i].tipo == 32:
                i_reagente_0 = self.specie.index(self.reazioni[i].reagenti[0])
                i_reagente_1 = self.specie.index(self.reazioni[i].reagenti[1])
                i_reagente_2 = self.specie.index(self.reazioni[i].reagenti[2])
                i_prodotto_0 = self.specie.index(self.reazioni[i].prodotti[0])
                i_prodotto_1 = self.specie.index(self.reazioni[i].prodotti[1])
                flusso = self.reazioni[i].costante * specie[i_reagente_0] * specie[i_reagente_1] * specie[i_reagente_2]
                if specie[i_reagente_0] <= 0:
                    flusso = 0
                if specie[i_reagente_1] <= 0:
                    flusso = 0
                if specie[i_reagente_2] <= 0:
                    flusso = 0

                delta[i_reagente_0] -= flusso
                delta[i_reagente_1] -= flusso
                delta[i_reagente_2] -= flusso
                delta[i_prodotto_0] += flusso
                delta[i_prodotto_1] += flusso
            elif self.reazioni[i].tipo == 210:
                i_reagente_0 = self.specie.index(self.reazioni[i].reagenti[0])
                i_prodotto_0 = self.specie.index(self.reazioni[i].prodotti[0])
                # qui ha senso un flusso negativo
                flusso = self.reazioni[i].costante * (specie[i_reagente_0] -
                                                      specie[i_prodotto_0]) / \
                                                      pow(self.volume, 1 / 3)
                delta[i_reagente_0] -= flusso
                delta[i_prodotto_0] += flusso
            else:
                print("ALLARME FN: ", self.reazioni[i].tipo)
                exit()

        for s in self.specie:
            if bool(s.inter["bufferizzata"]) is True:
                delta[self.specie.index(s.nome)] = 0

        # Controllo specie negativa ---> mando specie a 0
        for s in self.specie:
            if specie[self.specie.index(s.nome)] + delta[self.specie.index(s.nome)] < 0:
                delta[self.specie.index(s.nome)] = -specie[self.specie.index(s.nome)]

        delta[-1] = dC

        return delta

    def duplicate(self, lipidi):
        """Effettua la duplicazione basandosi sulla quantità di lipide passata"""
        qnt_dup = lipidi / 2
        self.volume = self.calc_volume(qnt_dup)
        return qnt_dup

    def fill_division_history(self, t, y):
        try:
            t_offset = list(self.division_history.keys())[-1]
        except IndexError:
            t_offset = 0
        self.division_history.update({t+t_offset: [s for s in y]})

    def fill_full_history(self, t, y):
        n_step = y.shape[1]
        try:
            t_offset = list(self.history.keys())[-1]
        except IndexError:
            t_offset = 0
        self.history.update({t[idx]+t_offset: [s[idx] for s in y] for idx in range(n_step)})

    def print_info(self):
        if len(self.history) != 0:
            self.print_full_history_to_file()
            self.print_history_graph()

        if len(self.division_history) != 0:
            self.print_division_file()
            self.print_division_graph()
        plt.show()

    def print_full_history_to_file(self):
        file = self.chem_file.split(".")[0] + "_out.txt"

        with open(file, "w") as f:
            f.write("T" + "\t" + str([s.nome for s in self.specie]).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "") + "\t" + "Container" + "\n")
            for t, l in self.history.items():
                f.write(str(t)+"\t"+str(l).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "")+"\n")

    def print_history_graph(self):
        x = list(self.history.keys())

        plt.figure(1)

        for s in self.specie:
            y = [l[self.specie.index(s)] for t, l in self.history.items()]
            plt.plot(x, y, label=s.nome)

        plt.xlabel('time (s)')

        plt.legend()
        # plt.show()

    def print_division_file(self):
        file = self.chem_file.split(".")[0] + "_division.txt"
        t_prev = 0
        with open(file, "w") as f:
            f.write("T" + "\t" + "T_prev_div" + "\t" + str([s.nome for s in self.specie]).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "") + "\t" + "Container" + "\n")

            for t, l in self.division_history.items():
                t_offset = t - t_prev
                f.write(str(t)+"\t"+str(t_offset)+"\t"+str(l).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "")+"\n")
                t_prev = t

    def print_division_graph(self):
        x = [i for i in range(len(self.division_history))]

        plt.figure(2)

        for s in self.specie:
            y = [l[self.specie.index(s)] for t, l in self.division_history.items()]
            plt.plot(x, y, label=s.nome)

        plt.xlabel('n. divisioni')

        plt.legend()
        # plt.show()

    def terminate(self, t, y):
        """Condizione di terminazione --> quando la quantità di lipide è raddoppiata"""
        if y[-1]/self.contenitore >= self.div_coeff:
            return 0
        return 1
    terminate.terminal = True

    def calc_volume(self, lipidi):
        return (1000 / 6) * pi * pow(self.membrane_thickness, 3) * \
               pow(pow((lipidi / (self.density * pi * pow(self.membrane_thickness, 3))) - 1 / 3, 1 / 2) - 1, 3)

    @staticmethod
    def read_conf(conf_file):
        with open(conf_file) as f:
            f.readline()    # riga contenente i titoli

            riga = f.readline()
            x0 = int(riga.split()[0])
            xn = int(riga.split()[1])
            n = int(riga.split()[2])
            n_div = int(riga.split()[3])
            div_coeff = int(riga.split()[4])
            print_division_file = int(riga.split()[5])
            print_history_file = int(riga.split()[6])

            return x0, xn, n, n_div, div_coeff, print_division_file, print_history_file

    @staticmethod
    def pr_sc_vett(sc, vett):
        return [specie(s.nome, sc * s.qnt, s.inter) for s in vett]

    @staticmethod
    def somma_vett(v1, v2):
        return [specie(v1[i].nome, v1[i].qnt + v2[i].qnt, v1[i].inter) for i in range(len(v1))]

    @staticmethod
    def copy(specie1):
        """ Deep copy """
        return [specie(s1.nome, s1.qnt, s1.inter) for s1 in specie1]

    def __str__(self):
        return f"Numero specie: {str(self.n_specie)}\nNumero reazioni: {self.n_reazioni}\nSpecie: {str(self.specie)}\n"
