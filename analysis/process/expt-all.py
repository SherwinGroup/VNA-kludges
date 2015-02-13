import numpy as np
import refltool as rt
from os import chdir 
from pylab import *


############
##import data

chdir("../data/expt2")
refl_expt2=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","half-sapph-no-ito_VNATrc.001","half-sapph-no-ito_VNATrc.002")


chdir("../expt3")
refl_expt3=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","half-sapph-yes-ito_VNATrc.001","half-sapph-yes-ito_VNATrc.002")


chdir("../expt4")
refl_expt4=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","full-coverage-ito-up_VNATrc.001","full-coverage-ito-up_VNATrc.002")


chdir("../expt5")
refl_expt5=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","full-coverage-ito-down_VNATrc.001","full-coverage-ito-down_VNATrc.002")


chdir("../expt6")
refl_expt6=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","si-ito-up_VNATrc.001","si-ito-up_VNATrc.002")


chdir("../expt7")
refl_expt7=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","si-ito-down_VNATrc.001","si-ito-down_VNATrc.002")


chdir("../sapph-old/clean")
refl_expt_o1=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","sapphire_VNATrc.001","sapphire_VNATrc.002")


chdir("../su8-bot")
refl_expt_o2=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","su8onsaph_VNATrc.001","su8onsaph_VNATrc.002")


chdir("../su8-top")
refl_expt_o3=rt.reflm_std("mirror_VNATrc.001","mirror_VNATrc.002","saff_VNATrc.001","saff_VNATrc.002")


chdir("../../../process")

#############################################


samps={}
freqbd={}
freqbd['wr3p4']=[228e9,311e9]
freqbd['wr2p2']=[330e9,490e9]
freqbd['wr1p5']=[490e9,680e9]


samps['half-sapph-no-ito']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.482e-3,window=100)
samps['half-sapph-no-ito'].add_trc('wr1p5', refl_expt2)


samps['half-sapph-yes-ito']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.482e-3,window=100)
samps['half-sapph-yes-ito'].add_trc('wr1p5', refl_expt3)


samps['full-coverage-ito-up']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.486e-3,window=100)
samps['full-coverage-ito-up'].add_trc('wr1p5', refl_expt4)


samps['full-coverage-ito-down']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.486e-3,window=100)
samps['full-coverage-ito-down'].add_trc('wr1p5', refl_expt5)


samps['si-ito-up']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.542e-3,window=100)
samps['si-ito-up'].add_trc('wr1p5', refl_expt6)


samps['si-ito-down']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.542e-3,window=100)
samps['si-ito-down'].add_trc('wr1p5', refl_expt7)


samps['sapph-clean']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.5e-3,window=100)
samps['sapph-clean'].add_trc('wr1p5', refl_expt_o1)


samps['sapph-su8-bot']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.5e-3,window=100)
samps['sapph-su8-bot'].add_trc('wr1p5', refl_expt_o2)


samps['sapph-su8-top']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.5e-3,window=100)
samps['sapph-su8-top'].add_trc('wr1p5', refl_expt_o3)


##############################################
##############################################

'''
plot(samps['sapph-clean'].trc_bd_fstack()[:,0], samps['sapph-clean'].efield_surf_abs())
plot(samps['sapph-clean'].trc_bd_fstack()[:,0], samps['sapph-su8-bot'].efield_surf_abs())
plot(samps['sapph-clean'].trc_bd_fstack()[:,0], samps['sapph-su8-top'].efield_surf_abs())
show()
'''

for samp in samps:
    np.savetxt(samp+"-efield.txt", np.hstack((
        samps[samp].trc_bd_fstack()[:,0],
        samps[samp].efield_surf_abs()
        ))
    )

