import matplotlib.pyplot as plt

def leggi_file_chimica():

    specie_nomi=[]
    q_specie=[]
    specie_inter=[]

    f=open("chem.txt","r")

    riga=f.readline()
    n_specie=eval(riga.split()[1])
    riga=f.readline()
    n_reazioni=eval(riga.split()[1])

    riga=f.readline()#riga vuota

    for i in range(n_specie):
        riga=f.readline()
        specie_nomi+=[riga.split()[0]]
        q_specie+=[eval(riga.split()[1])]
        specie_inter+=[{"buff":eval(riga.split()[2]),"cont":eval(riga.split()[3])}]     # 3°-4° colonna buff: bufferizzata  cont: contatto con contenitore

    #CONTROLLARE CHE NON CI SIANO NOMI RIPETUTI

    riga=f.readline()#riga vuota

    reazioni=[]
    for i in range(n_reazioni):
        reazione={}
        riga=f.readline()
        reazione["tipo"]=eval(riga.split()[0])
        if reazione["tipo"]==21:
            reazione["reagenti"]=[specie_nomi.index(riga.split()[1]),specie_nomi.index(riga.split()[3])]
            reazione["prodotti"]=[specie_nomi.index(riga.split()[5])]
            reazione["costante"]=eval(riga.split()[7])
        elif reazione["tipo"]==12:
            reazione["reagenti"]=[specie_nomi.index(riga.split()[1])]
            reazione["prodotti"]=[specie_nomi.index(riga.split()[3]),specie_nomi.index(riga.split()[5])]
            reazione["costante"]=eval(riga.split()[7])
        elif reazione["tipo"]==22:
            reazione["reagenti"]=[specie_nomi.index(riga.split()[1]),specie_nomi.index(riga.split()[3])]
            reazione["prodotti"]=[specie_nomi.index(riga.split()[5]),specie_nomi.index(riga.split()[7])]
            reazione["costante"]=eval(riga.split()[9])
        elif reazione["tipo"]==1:
            reazione["reagenti"]=[]
            reazione["prodotti"]=[specie_nomi.index(riga.split()[3])]
            reazione["costante"]=eval(riga.split()[5])
        else:
            print("ALLARMEleggi_file_chimica Tipo reazione: ",reazione["tipo"])

        reazioni+=[reazione]
        #print(riga.split())
        print(reazione)

    f.close()

    return (specie_nomi,q_specie,specie_inter,reazioni)


# function to be solved, ATM t è inutile
def fn(t,x,r):
    ''' Restituisce un vettore che indica le nuove
        quantità a seguito delle reazioni'''
    delta=[0]*len(x)
    for i in range(len(r)):
        if r[i]["tipo"]==12:
            flusso=r[i]["costante"]*x[r[i]["reagenti"][0]]
            delta[r[i]["reagenti"][0]]-=flusso
            delta[r[i]["prodotti"][0]]+=flusso
            delta[r[i]["prodotti"][1]]+=flusso
        elif r[i]["tipo"]==21:
            flusso=r[i]["costante"]*x[r[i]["reagenti"][0]]*x[r[i]["reagenti"][1]]
            delta[r[i]["reagenti"][0]]-=flusso
            delta[r[i]["reagenti"][1]]-=flusso
            delta[r[i]["prodotti"][0]]+=flusso
        elif r[i]["tipo"]==22:
            flusso=r[i]["costante"]*x[r[i]["reagenti"][0]]*x[r[i]["reagenti"][1]]
            delta[r[i]["reagenti"][0]]-=flusso
            delta[r[i]["reagenti"][1]]-=flusso
            delta[r[i]["prodotti"][0]]+=flusso
            delta[r[i]["prodotti"][1]]+=flusso
        elif r[i]["tipo"]==1:
            flusso=r[i]["costante"]
            delta[r[i]["prodotti"][0]]+=flusso
        else:
            print("ALLARME FN: ",r[i]["tipo"])
            exit()

    # RESITTUISCE UN VETTORE CON LE NUOVE QUANTITà
    return delta

"""
# function to be solved
def f(t,x):

    a=-Kad*x[0]*x[2]+Kai*x[3]+INA
    b=-Kbd*x[1]*x[3]+INB
    c=2*Kbd*x[1]*x[3]-Kad*x[0]*x[2]+Kai*x[3]
    compl=Kad*x[0]*x[2]-Kai*x[3]-Kbd*x[1]*x[3]
    return (a,b,c,compl)
"""

def pr_sc_vett(sc,vett):
    vett=list(vett)
    for i in range(len(vett)):
        vett[i]=sc*vett[i]
    return vett

def somma_vett(v1,v2):
    # vett = v1
    for i in range(len(v1)):
        v1[i]=v1[i]+v2[i]
    return v1

# RK-4 method
#   rk4(t0,x0,tn,step,reazioni,specie_inter,specie_nomi)
def rk4(x0,y0,xn,n,r,specie_inter,specie_nomi):

    # Calculating step size
    h = (xn-x0)/n

    time=[]
    ya=[]
    yb=[]
    yc=[]
    ycompl=[]
    time=time+[x0]
    ya=ya+[y0[0]]
    yb=yb+[y0[1]]
    yc=yc+[y0[2]]
    ycompl=ycompl+[y0[3]]

    y_iniziale=y0[:]# salvo le codizioni iniziali

    f=open("proto_pan1_out.txt","w")
    f.write("Time"+str(specie_nomi).replace(",","\t").replace("[","\t").replace("]","\n"))

    for i in range(n):
        # stampa ogni 100 passi, dovrà stampare ogni tot secondi
        if(i%100==0):
            print(y0)
            # x0 = tempo, y0 = lista quantità delle specie
            f.write(str(x0)+str(y0).replace(",","\t").replace("[","\t").replace("]","\n"))


        k1=pr_sc_vett(h,fn(x0, y0,r))
        k2=pr_sc_vett(h,fn(x0+h/2,somma_vett(y0,pr_sc_vett(0.5,k1)),r))
        k3=pr_sc_vett(h,fn(x0+h/2,somma_vett(y0,pr_sc_vett(0.5,k2)),r))
        k4=pr_sc_vett(h,fn(x0+h,somma_vett(y0,k3),r))

        k=somma_vett(k1,pr_sc_vett(2,k2))
        k=somma_vett(k,pr_sc_vett(2,k3))
        k=somma_vett(k,k4)
        k=pr_sc_vett(1/6,k)

        yn=somma_vett(y0,k)

        #print('%.4f\t%.4f\t%.4f'%(x0,y0[0],y0[1]) )

        #k1 = h * (f(x0, y0))
        #k2 = h * (f((x0+h/2), (y0+k1/2)))
        #k3 = h * (f((x0+h/2), (y0+k2/2)))
        #k4 = h * (f((x0+h), (y0+k3)))
        #k = (k1+2*k2+2*k3+k4)/6
        #yn = y0 + k
        #print('%.4f\t%.4f\t%.4f'% (x0,y0,yn) )
        #print('-------------------------')

        #rimetto a posto le specie bufferizzate
        for i in range(len(y_iniziale)):
            if specie_inter[i]["buff"]==1:
                #print("i: %d\ty0: %f\tyn: %f"%(i,y_iniziale[i],yn[i]))
                yn[i]=y_iniziale[i]


        y0 = yn
        x0 = x0+h

        time=time+[x0]
        ya=ya+[yn[0]]
        yb=yb+[yn[1]]
        yc=yc+[yn[2]]
        ycompl=ycompl+[yn[3]]

    f.write(str(x0)+str(y0).replace(",","\t").replace("[","\t").replace("]","\n"))
    f.close()

    print('Finale: %.4f\t%.4f\t%.4f\n'%(x0,y0[0],y0[1]) )

    f1 = plt.figure(1)
    #plt.subplot(211)

    plt.plot(time, ya, color = 'red', label='ya')
    plt.plot(time, yb, color = 'blue', label='yb')
    plt.plot(time, yc, marker = "o",color = 'green', label='yc')
    plt.plot(time, ycompl, color = 'orange', label='ycompl')

    plt.legend(numpoints=2)

    plt.show()

#(specie_nomi,x,specie_inter,reazioni) --> x è la lista contenente le quantità per ogni specie
(specie_nomi,q_specie,specie_inter,reazioni)=leggi_file_chimica()
"""
for i in range(len(x0)):
    print("%s\t%f\t%d\t%d\n"%(specie_nomi[i],x0[i],specie_inter[i]["buff"],specie_inter[i]["cont"]))
"""
"""
for i in range(len(reazioni)):
    print("R%d\t"%i,end="\t")
    print(reazioni[i])
"""

# !!! QUESTI PARAMETRI ANDRANNO PRESI DA UN FILE ESTERNO !!!

# t0 = float(input('t0 = '))
#
# print('Tempo finale: ')
# tn = float(input('tn = '))
#
# print('Enter number of steps:')
# step = int(input('Number of steps = '))

# RK4 method call
rk4(0,q_specie,20,10000,reazioni,specie_inter,specie_nomi)




# RK-4 method python program

# function to be solved
def f(x,y):
    return y
    #MVreturn x+y

# or
# f = lambda x: x+y
