import numpy as np
import refltool as rt
from os import chdir
from pylab import *

chdir("../data/expt1")

# structure of this data file is 
# (freq) (log mag ~ reflectivity) (Re(s11)) (Im(s11))
data=np.loadtxt("quartz-transm_VNATrc.002")

rv=data[:,2] #??
freq=data[:,0]

#rv= np.abs(data[:,3])**2 + np.abs(data[:,4])**2

plot(freq,rv)
#plot(freq,np.log(rv)/np.log(10)) 
show()

