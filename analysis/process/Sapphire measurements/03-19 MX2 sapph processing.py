# -*- coding: utf-8 -*-
"""
Created on Mon Feb 23 13:07:08 2015

@author: hbanks
"""
import os
import numpy as np
import matplotlib.pyplot as plt
import refltool as rt
#C:\Users\hbanks\Documents\GitHub\VNA-kludges\analysis\process\Sapphire measurements

substrate = rt.reflm_std("WR3o4 Sapphire MX2 substrate\mirrorr_VNATrc.001", 
                         "WR3o4 Sapphire MX2 substrate\mirrorr_VNATrc.002", 
                         "WR3o4 Sapphire MX2 substrate\sapph_VNATrc.001", 
                         "WR3o4 Sapphire MX2 substrate\sapph_VNATrc.002")

samps = {}
freqbd = {}
freqbd['wr3p4'] = [228e9,311e9]
freqbd['wr1p5'] = [490e9,680e9]

samps['substrate'] = rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.482e-3,window=100)
samps['substrate'].add_trc('wr3p4', substrate)


for elem in samps.values():
    a = rt.savitzky_golay(elem.efield_surf_abs(), 51, 2, deriv=0, rate=1)
    plt.plot(elem.trc_bd_fstack()[:,0]*1e-9, a)
plt.show()

for samp in samps:
    np.savetxt(samp+"-efield.txt", np.hstack((
        samps[samp].trc_bd_fstack()[:,0].reshape((len(samps[samp].trc_bd_fstack()[:,0]),1)),
        samps[samp].efield_surf_abs().reshape((len(samps[samp].trc_bd_fstack()[:,0]),1))
        ))
    )


