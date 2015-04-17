import matplotlib.pyplot as plt
from matplotlib.ticker import AutoMinorLocator
import numpy as np

x=np.linspace(100,200,num=1000)
y=85.5*np.sin(x/5.)/x*1e9
y2=x**2/20000.*1e9
y3=-y/y2*1e9
y4=1/y2*1e18
y5=np.cos(x/7.)*y2
y6=y5*np.sin(x/8.)

fig, ax = plt.subplots()
plt.minorticks_on()

plt.plot(x,y)
plt.plot(x,y2)
plt.plot(x,y3)
plt.plot(x,y4)
plt.plot(x,y5)
plt.plot(x,y6)




#labelpad
plt.xlabel("Time ($\sigma_{22}$)",labelpad=3)
plt.ylabel("WWW",labelpad=3)
#plt.yticks((-.951e9,0.,.951e9,1.902e9))
#plt.ylim(-.951e9,1.902e9)
#plt.yticks((-1e9,0.,1e9,2e9))
#plt.xticks((100,140,180),minor=True)

ax.xaxis.set_minor_locator(AutoMinorLocator(4))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))

plt.savefig("tesfig.pdf")