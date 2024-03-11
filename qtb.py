import random as rd 
import numpy as np 


from functions import *
n=int(pas/2-1)
def nb_aleatoire(nb):
    return [complex(rd.gauss(0,1)) for i in range(nb)]

    
def Theta(w,temp,type):
    beta=1/(kb*temp)
    #print(beta*hb*w)
    if w<6*10**(15):

        if type=="q":
            x=hb*w/2 +hb*w*(np.exp(beta*hb*w)-1)**(-1)
            return x
        else: return kb*temp
    else: 
        return 0
    


def IA(temp,n,type):
    dw=2*np.pi/pas/dt 
    IA=[0]+[2*m*gamma*Theta(i*dw,temp,type) for i in range(1,n+1)]
    return IA

def A_tilde(Ntilde,Mtilde,IA):
    At=np.sqrt(pas*dt/2)*np.array([np.sqrt(IA[i])*(Ntilde[i]+Mtilde[i]*1j) for i in range(len(Ntilde))])
    return At



def Atilde_sym(At):
    At=At.tolist()
    At.append(0)
    for i in range(1,n+1):
        At.append(np.conjugate(At[n+1-i]))
    return At




def A_fft(A):
    A_values=np.fft.ifft(np.array(A)/dt/pas)
    return np.real(A_values)*pas

# def real_tf(A):
#     x=0
#     B=[]
#     for i in range(pas):
#         x=0
#         for k in range(pas):
#         x+=A[k]*np.exp(-2*1j*np.pi*i*k/)


def inverse_fourier_transform(signal):
    N = len(signal)
    result = []

    for n in range(N):
        print(n)
        sum_real = 0
        sum_imaginary = 0

        for k in range(N):
            angle = 2 * np.pi * k * n / N
            sum_real += signal[k] * np.cos(angle)
            sum_imaginary += signal[k] * np.sin(angle)

        result.append((sum_real + sum_imaginary * 1j) / N)

    return result



def F_qtb(temp,type="q"):
    At=A_tilde(nb_aleatoire(n+1),nb_aleatoire(n+1),IA(temp,n,type))
    A_sym=Atilde_sym(At)
    return A_fft(A_sym)

def F_qtb_exacte(temp,type="q"):
    At=A_tilde(nb_aleatoire(n+1),nb_aleatoire(n+1),IA(temp,n,type))
    A_sym=Atilde_sym(At)
    return inverse_fourier_transform(A_sym)


# At=A_tilde(nb_aleatoire(n+1),nb_aleatoire(n+1),IA(T,n,"c"))
# A_sym=Atilde_sym(At)
# print(A_sym)