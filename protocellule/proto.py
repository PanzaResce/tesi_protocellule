import os
import re
import warnings
import matplotlib.pyplot as plt
from scipy.integrate import solve_ivp
from scipy.integrate._ivp.ivp import METHODS
from math import pi
from protocellule import scipy_patch
from protocellule.specie import specie
from protocellule.reazione import reazione


class Proto:
    """ Per come è concepita la classe proto, tutte
        le informazioni che vengono memorizzate rappresentanto
        lo stato iniziale della cellula """

    def __init__(self, conf_file="conf.txt", chem_file="chem.txt"):

        self.specie = list()
        self.reazioni = list()
        # self.history = dict()
        self.history = list()
        self.division_history = dict()
        # self.flow_history = list()

        # Event-related attribute
        self.t_abs = 0
        self.t_prev = 0
        self.t_hist = list()
        self.events = list()

        self.conf_file = conf_file
        self.chem_file = chem_file
        self.output_dir = os.path.join("./", chem_file.split(".")[0] + "_output")

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

                # this is for optimization purpose, the part (36pi^1/3 / membrane) is reaction_constant through all the program
                if tipo == "210":
                    vett_reazione[2] = float(vett_reazione[2]) * pow(36 * pi, 1 / 3) / self.membrane_thickness

                self.reazioni.append(reazione(tipo, vett_reazione))

            # Lettura eventi
            res = input("La simulazione prevede degli eventi? (y/n)\n")
            if res == "y":
                print("Gli eventi saranno presi dalla directory ./eventi\n")
                eventi = list()
                for filename in os.scandir("./eventi"):
                    if filename.is_file():
                        # print(filename.path)
                        with open(filename) as f:
                            riga = f.readline()
                            t = float(riga.split()[1])
                            reaction_const = dict()
                            species_concentration = dict()
                            f.readline()   #riga vuota
                            f.readline()   #riga intestazioni

                            for line in f:
                                # if it's a concentration change event
                                if line.split()[0] in self.specie:
                                    specie_tipo = line.split()[0]
                                    concentrazione = line.split()[1]
                                    species_concentration[specie_tipo] = float(concentrazione)
                                    continue
                                # if it's a reaction constant change event
                                n_reazione = int(line.split()[0])
                                if n_reazione > len(self.reazioni):
                                    warnings.warn(f"Numero di reazione non valido, le reazioni sono {len(self.reazioni)}", UserWarning)
                                    self.events = list()
                                    return
                                reaction_const[n_reazione] = float(line.split()[1])

                            evento = (t, reaction_const, species_concentration)
                            eventi.append(evento)
                # files are not opened tidily, so the list needs to be ordered by the event time
                self.events = Proto.sort_list(eventi)
                print(self.events)

        if not os.path.exists(self.output_dir):
            os.mkdir(self.output_dir)

    def simula(self):
        (x0, xn, n, n_div, div_coeff, print_div, print_hist, integration_method, lower_threshold) = Proto.read_conf(self.conf_file)
        self.div_coeff = div_coeff

        t_span = [x0, xn]

        y0 = [s.qnt for s in self.specie]

        # Nel vettore passato all'integratore gli ultimi due elementi sono rispettivamente
        # il lipide vecchio e il lipide nuovo, all'inizio sono uguali
        # Servono per calcolare il volume vecchio e il volume nuovo
        y0.append(self.contenitore)
        y0.append(self.contenitore)

        for i in range(n_div):
            sol1 = solve_ivp(self.fn, t_span, y0, method=integration_method, events=[self.terminate, self.check_event])

            # va fatto un ulteriore controllo sulle quantità negative
            for step in sol1.y:
                for idx in range(len(step)):
                    if step[idx] < lower_threshold:
                        step[idx] = lower_threshold

            n_step = sol1.y.shape[1]
            out = [s[n_step - 1] for s in sol1.y]

            print(f"N. divisione: {i+1}, {out[-1]}, {sol1.t[-1]}, {self.t_abs}")

            if print_div:
                self.fill_division_history(sol1.t[n_step-1], out)
            if print_hist:
                self.fill_full_history(sol1.t, sol1.y, gen=i)

            new_qnt = self.duplicate(out[-1])

            out[-1] = new_qnt
            out[-2] = new_qnt

            # apply events here
            for ev in self.events:
                if ev[0] <= self.t_abs:
                    self.apply_event(ev, out)     # side-effect on out
                    print(f"TRIGGER: {len(self.events)-1} events remaining")

            # delete happened event
            self.events = [ev for ev in self.events if ev[0] > self.t_abs]

            y0 = out

        print("END")

        #Test
        # c = 0
        # for i in range(len(self.t_hist)):
        #     if i != 0:
        #         if self.t_hist[i-1] > self.t_hist[i]:
        #             print(f"Errore {i}")
        #             c+=1
        # print(f"Errori: {c}")

        #Test
        print(f"{len(self.events)} eventi non avvenuti")

    def fn(self, arr_t, specie):
        """ Restituisce un vettore che indica la variazione delle nuove concentrazioni a seguito delle reazioni
         Parametri:
             t: lista
                t[0] : istante di tempo attuale
                t[1] : dt
             specie : lista che contiene le concentrazioni delle specie
         """

        # if t > self.dt[1]:
        #     dt = t - self.dt[1]
        # elif t == self.dt[1]:
        # elif t < self.dt[1]:

        # Unpack arr_t
        try:
            t = arr_t[0]
            dt = arr_t[1]
        except (TypeError, IndexError):
            t = arr_t
            dt = 1

        # Calcolo volume attuale con quantità di lipide attuale
        self.volume = self.calc_volume(specie[-1])
        # Calcolo volume vecchio con quantità di lipide precedente
        Vv = self.calc_volume(specie[-2])

        delta = [0] * len(specie)

        if t != 0:
            v_rapp = Vv / self.volume

            # ricalcolo concentrazioni sostanze con rapporto (volume vecchio / volume nuovo)
            for idx, s in enumerate(specie):
                # if idx != len(specie)-1 and bool(self.specie[idx].inter["bufferizzata"]) is not True:
                if idx != len(specie)-1 and idx != len(specie)-2 and bool(self.specie[idx].inter["bufferizzata"]) is not True:
                    delta[idx] -= specie[idx]*(1-v_rapp) / dt

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

        # calcolo variazione quantità di lipide con nuovo volume
        dC = sum([s.inter["boundary"] * specie[self.specie.index(s)] for s in self.specie]) * Vv

        delta[-1] = dC
        delta[-2] = (specie[-1] - specie[-2]) / dt

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

    def fill_full_history(self, t, y, gen):
        n_step = y.shape[1]
        try:
            # t_offset = list(self.history.keys())[-1]
            t_offset = self.history[-1][0]
        except IndexError:
            t_offset = 0
        # self.history.update({t[idx]+t_offset: [s[idx] for s in y] for idx in range(n_step)})
        local_history = list()
        for idx in range(n_step):
            local_history.append([t[idx] + t_offset] + [s[idx] for s in y])

        self.history += local_history

        file = self.chem_file.split(".")[0] + f"_gen{gen}_out.part.txt"
        self.print_full_history_to_file(file, local_history)

    def print_final_info(self):
        if len(self.history) != 0:
            file = self.chem_file.split(".")[0] + "_out.txt"
            self.print_full_history_to_file(file, self.history)
            self.print_history_graph()

        if len(self.division_history) != 0:
            self.print_division_file()
            self.print_division_graph()
        plt.show()

        # clear intermediate generation's file
        for file in os.listdir(self.output_dir):
            if re.search(r".*gen\d+_out\.part\.txt$", file):
                os.remove(os.path.join(self.output_dir, file))

    def print_full_history_to_file(self, filename, dict):

        with open(os.path.join(self.output_dir, filename), "w") as f:
            f.write("T" + "\t" + str([s.nome for s in self.specie]).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "") + "\t" + "Container" + "\n")
            # for t, l in self.history.items():
            #     f.write(str(t)+"\t"+str(l).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "")+"\n")
            for step in dict:
                f.write(str(step).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "") + "\n")

    def print_history_graph(self):
        # x = list(self.history.keys())
        x = [step[0] for step in self.history]

        plt.figure(1)

        for s in self.specie:
            # y = [l[self.specie.index(s)] for t, l in self.history.items()]
            y = [step[self.specie.index(s)+1] for step in self.history]
            plt.plot(x, y, label=s.nome)

        plt.xlabel('time (s)')

        plt.legend()
        plt.yscale('log')
        plt.ylim(bottom=1e-10)
        # plt.show()

    def print_division_file(self):
        file = self.chem_file.split(".")[0] + "_division.txt"
        t_prev = 0
        with open(os.path.join(self.output_dir, file), "w") as f:
            f.write("T" + "\t" + "T_prev_div" + "\t" + str([s.nome for s in self.specie]).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "") + "\t" + "Container" + "\n")

            for t, l in self.division_history.items():
                t_offset = t - t_prev
                f.write(str(t)+"\t"+str(t_offset)+"\t"+str(l).replace(",", "\t").replace("[", "").replace("]", "").replace("'", "")+"\n")
                t_prev = t

        # with open("prova.txt", "w") as f:
        #     for t in self.event_history:
        #         f.write(str(t)+"\n")

    def print_division_graph(self):
        x = [i for i in range(len(self.division_history))]

        plt.figure(2)

        for s in self.specie:
            y = [l[self.specie.index(s)] for t, l in self.division_history.items()]
            plt.plot(x, y, label=s.nome)

        plt.xlabel('n. divisioni')

        plt.legend()
        plt.yscale('log')
        plt.ylim(bottom=1e-10)
        # plt.show()

    def terminate(self, t, y):
        """
        Condizione di terminazione --> quando la quantità di lipide supera una soglia prestabilita
        Questa soglia dipende dal coefficiente passato tramite file conf
        Ex: div_coeff = 2 --> termina quando quantità lipide raddoppia rispetto a quella iniziale
        """
        if y[-1]/self.contenitore >= self.div_coeff:
            return 0
        return 1
    terminate.terminal = True

    def check_event(self, t, y):
        if t == 0:
            self.t_prev = 0

        prev_t = float("inf")

        # Calcolo tempo assoluto, ovveri i tempi in cui l'integratore chiama l'evento
        if t > self.t_prev:
            self.t_abs = self.t_abs + (t - self.t_prev)
            self.t_prev = t
            try:
                prev_t = self.t_hist[-1]
            except IndexError:
                prev_t = float("inf")
            self.t_hist.append(self.t_abs)

        # if len(self.events) <= 0:
        #     return 0

        # il primo elemento è il prossimo evento che deve avvenire
        # t_evento = self.events[0][0]
        #
        # if prev_t <= t_evento < self.t_abs:
        #     evento = self.events.pop(0)
        #     self.apply_event(evento, y)     # side-effect on y
        #     print("TRIGGER")

        return 0
    check_event.terminal = False

    # works by side-effect (specially on y)
    def apply_event(self, event, y):
        for n_reazione, const in event[1].items():
            if str(self.reazioni[n_reazione].tipo) == "210":
                const = const * pow(36 * pi, 1 / 3) / self.membrane_thickness
            self.reazioni[n_reazione].costante = const

        for s, concentration in event[2].items():
            y[self.specie.index(s)] = concentration

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
            integration_method = riga.split()[7]
            if integration_method not in METHODS.keys():
                warnings.warn(f"Method {integration_method} not valid, switching by default to LSODA method", UserWarning)
                integration_method = 'LSODA'
            lower_threshold = float(riga.split()[8])

            return x0, xn, n, n_div, div_coeff, print_division_file, print_history_file, integration_method, lower_threshold

    @staticmethod
    def sort_list(list):
        for i in range(1, len(list)):
            j = i
            while j > 0 and list[j-1][0] > list[j][0]:
                app = list[j]
                list[j] = list[j-1]
                list[j-1] = app
                j = j-1

        return list

    def __str__(self):
        return f"Numero specie: {str(self.n_specie)}\nNumero reazioni: {self.n_reazioni}\nSpecie: {str(self.specie)}\n"
