import random as rd
import numpy as np
import matplotlib.pyplot as plt
import re

class Atome:

    def __init__(self):
            self.x = [a,a,a]
            self.v = [1,0,0]




P,T,pas=0,0,0
with open('assets/data/variables.txt', 'r') as file:
        pattern = r'(\w+)=(\d+)'
        # Lire le contenu du fichier ligne par ligne
        for line in file:
            # Utiliser la regex pour chercher des correspondances dans chaque ligne
            matches = re.findall(pattern, line)
            # Si des correspondances sont trouvées, les assigner comme variables globales
            if matches:
                for match in matches:
                    variable, valeur = match
                    # Utilisation de la fonction globals() pour assigner comme variable globale
                    globals()[variable] = int(valeur)
m=1.6726*10**(-27)
ev=1.6*10**(-19)
V0=0.3*ev
kb=1.380649 *10**(-23)
hb=6.626*10**(-34)/2
hb/=np.pi
a=0.3*10**(-10)
e=a/100
w0 = np.sqrt((8*V0/(a**2*m)))

wP=P*kb*T/hb
dt=1/1000*2*np.pi/w0 
dt=10**(-16)  #attention 17 nrmlt
gamma=1/(1000*dt)
i=0
K=P**2*m*kb**2*T**2/hb**2/100
atome={}
def Vp(x):
    return V0*(1-(x/a)**2)**2

def Vpprime (x) :
    V_prime = (4*V0/a**2) * ((x**3/a**2) - x)
    return (V_prime)
def Vsin(x):
    return V0**2*np.sin(2*np.pi*x/a)**2

def Vsinprime(x):
    return 4*np.pi/a*V0**2*np.sin(2*np.pi*x/a)*np.cos(2*np.pi*x/a)



def init(N):   #creates N atoms with the same caracteristics
    for i in range(N):
        atome[i]=Atome()

def voisins(i,N): # gives the neighbours of a givin atom
    if N==1:
        return (0,0)
    if i==N-1:
        return (N-2,0)
    elif i==0:
        return (N-1,1)
    return i-1,i+1



def main_RPMD(etape,N,V,Vprime):  #integration of the motion equation, 
    for i in range(N):
        F = - 1/N*Vprime(atome[i].x[-1]) - m*gamma*atome[i].v[-1] + np.sqrt ((2*m*gamma*kb*T)/dt)*rd.gauss(0,1)
        pp,pn=voisins(i,N)
        Fk=-K*(2*atome[i].x[2+etape]-atome[pp].x[2+etape]-atome[pn].x[2+etape])
        F+=Fk
        pos = 2*atome[i].x[-1] - atome[i].x[-2]+dt**2 * F/m
        vit = (3*atome[i].x[-1]-4*atome[i].x[-2]+atome[i].x[-3])/(2*dt)
        atome[i].x.append(pos)
        atome[i].v.append(vit)

def simulation_RPMD(N,V,Vprime): #main function which computes the positions of the P particles in a given distribution of potential
    init(N)
    for i in range(pas):
        if i%100000==0:
            print(i)
        main_RPMD(i,N,V,Vprime)

def main_QTB(etape,N,V,Vprime,A):  #integration of the motion equation, 
    for i in range(N):
        F = - 1/N*Vprime(atome[i].x[-1]) - m*gamma*atome[i].v[-1] + A[etape]
        pos = 2*atome[i].x[-1] - atome[i].x[-2]+dt**2 * F/m
        vit = (3*atome[i].x[-1]-4*atome[i].x[-2]+atome[i].x[-3])/(2*dt)
        atome[i].x.append(pos)
        atome[i].v.append(vit)

def simulation_QTB(N,V,Vprime,A): #main function which computes the positions of the P particles in a given distribution of potential
    init(N)
    print(max(A))
    for i in range(pas):
        if i%100000==0:
            print(i)
        main_QTB(i,N,V,Vprime,A)


def Veff(B,V):
    x=0
    for i in B:
        x+=V(i)/kb/T
    return x/len(B)

def Vmoy_for_P(N,V):
    Vmoy=[0]
    for i in range(pas):
        Vmoy.append((Vmoy[-1]*(i+1)+Veff([atome[j].x[-1] for j in range(N)],V))/(i+2))
    return Vmoy[-1]


def vitesse_moy2(N):
    v2=0
    for i in range(N):
        for j in range(pas):
            v2+=(atome[i].v[j])**2

    return v2/N/pas

'''
simulation_RPMD(P,Vp,Vpprime)
print(vitesse_moy2(P)*m/kb)
'''