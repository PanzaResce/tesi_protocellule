import time
from protocellule.proto import proto

if __name__ == "__main__":

    file = input("Inserisci nome file(file.txt): ")

    start_time = time.time()

    p = proto("conf.txt", file+".txt")
    # p = proto()
    p.simula()

    print("--- %s seconds ---" % (time.time() - start_time))

    p.print_info()

    # visual.draw(p)
