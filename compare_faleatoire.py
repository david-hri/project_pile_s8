from functions import *
from qtb import *



R1=[np.sqrt ((2*m*gamma*kb*T)/dt)*rd.gauss(0,1)]
for i in range(pas-1):
    R1.append(np.sqrt ((2*m*gamma*kb*T)/dt)*rd.gauss(0,1))


R2=F_qtb(T,"q")


plt.plot([i for i in range(500)],R1[:500])
plt.plot([i for i in range(500)],R2[:500])
plt.show()
