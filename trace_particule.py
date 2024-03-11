
from functions import *
from qtb import *

print(0.073*ev/(0.5*kb))
print(dt)
print(kb)

simulation_QTB(1,Vp,Vpprime,F_qtb(1000,"q"))
plt.plot([i for i in range(len(atome[0].x))],atome[0].x)

plt.show()

print((vitesse_moy2(1)*m/kb))

'''
simulation_QTB(1,Vp,Vpprime,A_fft())

simulation_RPMD(P,Vp,Vpprime)
plt.plot([i for i in range(len(atome[0].x))],atome[0].x)
print(Vmoy_for_P(P,Vp))




 

'''