import time
from protocellule.proto import Proto

if __name__ == "__main__":

    file = input("Inserisci nome file(file.txt): ")

    start_time = time.time()

    p = Proto("conf.txt", file+".txt")
    # p = proto()
    p.simula()

    print("--- %s seconds ---" % (time.time() - start_time))

    p.print_final_info()

    # visual.draw(p)
