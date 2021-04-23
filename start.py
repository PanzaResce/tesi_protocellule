from protocellule.proto import proto

if __name__ == "__main__":
    file = input("Inserisci nome file(file.txt): ")
    p = proto("conf.txt", file+".txt")
    # p = proto()
    p.simula()
    p.print_to_file()
    p.print_division_file()
    p.print_graph()
