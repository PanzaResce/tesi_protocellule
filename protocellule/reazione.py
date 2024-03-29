class Reazione:

    def __init__(self, tipo, reazione):
                
        self.tipo = int(tipo)

        self.reagenti = list()
        self.prodotti = list()
        self.costante = float()

        if self.tipo == 21:     # creazione complesso
            self.reagenti = [reazione[0], reazione[1]]
            self.prodotti = [reazione[2]]
            self.costante = float(reazione[3])
        elif self.tipo == 12:   # scissione complesso
            self.reagenti = [reazione[0]]
            self.prodotti = [reazione[1], reazione[2]]
            self.costante = float(reazione[3])
        elif self.tipo == 22:   # scambio (complesso in prodotto)
            self.reagenti = [reazione[0], reazione[1]]
            self.prodotti = [reazione[2], reazione[3]]
            self.costante = float(reazione[4])
        elif self.tipo == 10:   # scomparsa
            self.reagenti = [reazione[0]]
            self.prodotti = []
            self.costante = float(reazione[2])
        elif self.tipo == 1:    # aggiunta
            self.reagenti = []
            self.prodotti = [reazione[1]]
            self.costante = float(reazione[2])
        elif self.tipo == 23:   # cleavage diretto
            self.reagenti = [reazione[0], reazione[1]]
            self.prodotti = [reazione[2], reazione[3], reazione[4]]
            self.costante = float(reazione[5])
        elif self.tipo == 32:   # condensazione
            self.reagenti = [reazione[0], reazione[1], reazione[2]]
            self.prodotti = [reazione[3], reazione[4]]
            self.costante = float(reazione[5])
        elif self.tipo == 210:  # passaggio di membrana
            self.reagenti = [reazione[0]]
            self.prodotti = [reazione[1]]
            self.costante = float(reazione[2])
