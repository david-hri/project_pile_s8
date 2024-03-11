from qtb import *


A=F_qtb_exacte(T,"q")
B=F_qtb(T,"q")


plt.plot([i for  i in range(200)],A[:200])
plt.plot([i for  i in range(200)],B[:200])
plt.show()