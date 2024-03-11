from functions import *
from qtb import *
import numpy as np

Nb_it=20
Tmin=10
Tmax=40000
T_table=np.logspace(np.log10(Tmin),np.log10(Tmax),Nb_it)



def calc_Tmd(temp):
    print(temp)
    simulation_QTB(1,Vp,Vpprime,F_qtb(temp,"q"))
    temp_md=(vitesse_moy2(1)*m/kb)
    print(temp_md)
    return temp_md

def tab_Tmd():
    Tmd=[]
    i=0
    for temp in T_table:
        i+=1
        print(i)
        Tmd.append(calc_Tmd(temp))
    return Tmd

Tmd=tab_Tmd()
plt.plot(T_table,Tmd)
plt.plot(T_table,T_table)
plt.show()

