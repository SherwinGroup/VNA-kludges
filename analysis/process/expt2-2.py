import numpy as np
import refltool as rt
from os import chdir 
from pylab import *


############
##import current data

#this one is half-sapph-no-ito
chdir("../data/expt2")
data_new=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","half-sapph-no-ito_VNATrc.001","half-sapph-no-ito_VNATrc.002")

###########
##import old data

chdir("../sapph-old")
data_old=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","sapphire_VNATrc.001","sapphire_VNATrc.002")



########
#inputs
#######

freq=data_new[:,0]

plot(freq,data_old[:,2:3])
plot(freq,data_new[:,2:3])
#plot(freq,abs(rt.f2c(data_new[:,2:4]))**2)
#plot(freq,abs(rt.f2c(data_old[:,2:4]))**2)
ylim(-1.5,1.5)
show()


'''
#incident e-field phase at surface
ph_i=0.
#..gives incident e-field, normalized to E_max = 1
e_inc=np.exp(1j*ph_i)


#reflected e-field
e_refl= e_inc*rt.f2c(data[:,2:4])[:,0]

#total field at surface is real part of ( incident e-field at surface plus reflected e-field at surface)


e_tot=e_inc+e_refl

plot(freq, e_tot.real)
plot(freq, e_refl.real)
ylim(-1.5,1.5)
show()
'''
