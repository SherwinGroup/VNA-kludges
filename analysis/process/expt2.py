import numpy as np
import refltool as rt
from os import chdir 
from pylab import *

#this one is half-sapph-no-ito

chdir("../data/expt2")

data=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","half-sapph-no-ito_VNATrc.001","half-sapph-no-ito_VNATrc.002")

########
#inputs
#######

freq=data[:,0]

#incident e-field phase at surface
ph_i=np.pi/2.
#..gives incident e-field, normalized to E_max = 1
e_inc=np.exp(1j*ph_i)


#reflected e-field
e_refl= e_inc*rt.f2c(data[:,2:4])[:,0]

#total field at surface is real part of ( incident e-field at surface plus reflected e-field at surface)

e_tot=e_inc+e_refl

plot(freq,e_tot.real)
plot(freq,e_tot.imag)
show()
