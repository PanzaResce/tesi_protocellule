import time
from protocellule.proto import proto

if __name__ == "__main__":

    file = input("Inserisci nome file(file.txt): ")

    start_time = time.time()

    p = proto("conf.txt", file+".txt")
    # p = proto()
    p.simula()

    print("--- %s seconds ---" % (time.time() - start_time))

    p.print_to_file()
    p.print_graph()

    p.print_division_file()

    # p.print_flow_file()

    # visual.draw(p)
