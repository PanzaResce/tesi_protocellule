import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from math import pi
from sample.specie import specie
from sample.reazione import reazione


class proto:
    """ Per come è concepita al momento la classe proto, tutte
        le informazioni che vengono memorizzate rappresentanto
        lo stato iniziale della cellula """

    def __init__(self, conf_file="conf.txt", chem_file="chem.txt"):

        self.specie = list()
        self.reazioni = list()
        self.history = list()

        self.conf_file = conf_file
        self.chem_file = chem_file

        with open(chem_file) as f:
            #Lettura proprietà cellula
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

            # Lettura specie --> nome, quantità, ...
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

                self.reazioni.append(reazione(tipo, vett_reazione))

    def simula(self):
        (x0, xn, n) = proto.read_conf(self.conf_file)

        t_span = [x0, xn]

        y0 = [s.qnt for s in self.specie]
        y0.append(self.contenitore)

        sol1 = solve_ivp(self.fn, t_span, y0, events=self.terminate)

        n_step = sol1.y.shape[1]
        out = [s[n_step - 1] for s in sol1.y]
        #self.rk4()

        print("---FINALE---")
        print(f"Vettore finale: {out}")
        print(f"Tempo: {sol1.t[-1]}")

        self.history = {sol1.t[i]: [s[i] for s in sol1.y] for i in range(n_step)}

        # x = range(n_step)
        # plt.plot(x, [x for x in self.history[0]], label="A")
        # plt.plot(x, [x for x in self.history[1]], label="B")
        # plt.plot(x, [x for x in self.history[2]], label="C")
        # plt.plot(x, [x for x in self.history[3]], label="*Compl")
        #
        # plt.legend()
        # plt.show()

        print(n_step)

    def fn(self, t, specie):
        """ Restituisce un vettore che indica le nuove
            quantità a seguito delle reazioni """

        # calcolo variazione quantità di lipide con vecchio volume
        dC = sum([s.inter["boundary"] * specie[self.specie.index(s)] for s in self.specie]) * self.volume

        if t != 0:
            # calcolo nuovo volume con nuova quantità
            new_volume = (1000/6)*pi*pow(self.membrane_thickness, 3)*pow(pow((specie[-1]/(self.density*pi*pow(self.membrane_thickness, 3)))-1/3, 1/2)-1, 3)
            v_rapp = self.volume / new_volume

            # aggiusto concentrazioni sostanze con rapporto volume vecchio / volume nuovo
            for idx, s in enumerate(specie):
                if idx != len(specie)-1 and bool(self.specie[idx].inter["bufferizzata"]) is not True:
                    specie[idx] *= v_rapp

            self.volume = new_volume

        delta = [0] * len(specie)

        # applico le reazioni (basandomi su vecchio volume)
        for i in range(self.n_reazioni):
            if self.reazioni[i].tipo == 12:
                flusso = self.reazioni[i].costante * specie[self.specie.index(self.reazioni[i].reagenti[0])]
                delta[self.specie.index(self.reazioni[i].reagenti[0])] -= flusso
                delta[self.specie.index(self.reazioni[i].prodotti[0])] += flusso
                delta[self.specie.index(self.reazioni[i].prodotti[1])] += flusso
            elif self.reazioni[i].tipo == 21:
                flusso = self.reazioni[i].costante * specie[self.specie.index(self.reazioni[i].reagenti[0])] * specie[self.specie.index(self.reazioni[i].reagenti[1])]
                delta[self.specie.index(self.reazioni[i].reagenti[0])] -= flusso
                delta[self.specie.index(self.reazioni[i].reagenti[1])] -= flusso
                delta[self.specie.index(self.reazioni[i].prodotti[0])] += flusso
            elif self.reazioni[i].tipo == 22:
                flusso = self.reazioni[i].costante * specie[self.specie.index(self.reazioni[i].reagenti[0])] * specie[self.specie.index(self.reazioni[i].reagenti[1])]
                delta[self.specie.index(self.reazioni[i].reagenti[0])] -= flusso
                delta[self.specie.index(self.reazioni[i].reagenti[1])] -= flusso
                delta[self.specie.index(self.reazioni[i].prodotti[0])] += flusso
                delta[self.specie.index(self.reazioni[i].prodotti[1])] += flusso
            elif self.reazioni[i].tipo == 10:
                flusso = self.reazioni[i].costante
                delta[self.specie.index(self.reazioni[i].reagenti[0])] -= flusso
            elif self.reazioni[i].tipo == 1:
                flusso = self.reazioni[i].costante
                delta[self.specie.index(self.reazioni[i].prodotti[0])] += flusso
            elif self.reazioni[i].tipo == 23:
                flusso = self.reazioni[i].costante * specie[self.specie.index(self.reazioni[i].reagenti[0])] * specie[self.specie.index(self.reazioni[i].reagenti[1])]
                delta[self.specie.index(self.reazioni[i].reagenti[0])] -= flusso
                delta[self.specie.index(self.reazioni[i].reagenti[1])] -= flusso
                delta[self.specie.index(self.reazioni[i].prodotti[0])] += flusso
                delta[self.specie.index(self.reazioni[i].prodotti[1])] += flusso
                delta[self.specie.index(self.reazioni[i].prodotti[2])] += flusso
            elif self.reazioni[i].tipo == 32:
                flusso = self.reazioni[i].costante * specie[self.specie.index(self.reazioni[i].reagenti[0])] * specie[self.specie.index(self.reazioni[i].reagenti[1])] * specie[self.specie.index(self.reazioni[i].reagenti[2])]
                delta[self.specie.index(self.reazioni[i].reagenti[0])] -= flusso
                delta[self.specie.index(self.reazioni[i].reagenti[1])] -= flusso
                delta[self.specie.index(self.reazioni[i].reagenti[2])] -= flusso
                delta[self.specie.index(self.reazioni[i].prodotti[0])] += flusso
                delta[self.specie.index(self.reazioni[i].prodotti[1])] += flusso
            elif self.reazioni[i].tipo == 210:
                flusso = ((pow(36*pi, 1/3) * self.reazioni[i].costante) / self.membrane_thickness) * (specie[self.specie.index(self.reazioni[i].reagenti[0])] - specie[self.specie.index(self.reazioni[i].prodotti[0])]) / pow(self.volume, 1/3)
                delta[self.specie.index(self.reazioni[i].reagenti[0])] -= flusso
                delta[self.specie.index(self.reazioni[i].prodotti[0])] += flusso
            else:
                print("ALLARME FN: ", self.reazioni[i].tipo)
                exit()

        for s in self.specie:
            if bool(s.inter["bufferizzata"]) is True:
                delta[self.specie.index(s.nome)] = 0

        delta[-1] = dC

        return delta

    def fill_history(self, vett):
        self.history.append(proto.copy(vett))

    def print_to_file(self, sec=None):
        file = self.chem_file.split(".")[0] + "_out.txt"
        with open(file, "w") as f:
            l = [el.nome for el in self.specie]
            f.write("Time\t"+str(l).replace(",", "\t").replace("[","").replace("]","").replace("'", "")+"\tLipide"+"\n")

            for t, l in self.history.items():
                f.write(str(t)+"\t"+str(l).replace(",", "\t").replace("[","").replace("]","").replace("'", "")+"\n")

    def print_graph(self, sec=None):

        x = list(self.history.keys())

        for s in self.specie:
            y = [l[self.specie.index(s)] for t, l in self.history.items()]
            plt.plot(x, y, label=s.nome)

        plt.legend()
        plt.show()

    def terminate(self, t, y):
        # termina quando [C] * volume_proto > Cmax
        # print(self.conf_file)
        # print(t)
        if y[-1]/self.contenitore >= 2.0:
            print(f"tempo finale: {t}")
            print(f"Vetttore finele: {y}")
            return 0
        return 1
    terminate.terminal = True

    @staticmethod
    def read_conf(conf_file):
        with open(conf_file) as f:
            f.readline()    # riga vuota

            riga = f.readline()
            x0 = int(riga.split()[0])
            xn = int(riga.split()[1])
            n = int(riga.split()[2])

            return x0, xn, n

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
