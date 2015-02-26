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

ITO_150_up = rt.reflm_std("WR3o4 150 sapph and blank sapph\ITOup_150nm_sapph_mirror_VNATrc.001", 
                          "WR3o4 150 sapph and blank sapph\ITOup_150nm_sapph_mirror_VNATrc.002", 
                          "WR3o4 150 sapph and blank sapph\ITOup_150nm_sapph_VNATrc.001", 
                          "WR3o4 150 sapph and blank sapph\ITOup_150nm_sapph_VNATrc.002")

ITO_150_down = rt.reflm_std("WR3o4 150 sapph and blank sapph\ITOdown_150nm_sapph_mirror_VNATrc.001", 
                            "WR3o4 150 sapph and blank sapph\ITOdown_150nm_sapph_mirror_VNATrc.002", 
                            "WR3o4 150 sapph and blank sapph\ITOdown_150nm_sapph_VNATrc.001", 
                            "WR3o4 150 sapph and blank sapph\ITOdown_150nm_sapph_VNATrc.002")

blank_3o4 = rt.reflm_std("WR3o4 150 sapph and blank sapph\\blank_sapph_mirror_VNATrc.001", 
                         "WR3o4 150 sapph and blank sapph\\blank_sapph_mirror_VNATrc.002", 
                         "WR3o4 150 sapph and blank sapph\\blank_sapph_VNATrc.001", 
                         "WR3o4 150 sapph and blank sapph\\blank_sapph_VNATrc.002")

blank_1o5 = rt.reflm_std("WR1o5 blank sapph\\blank_sapph_mirror_VNATrc.001", 
                         "WR1o5 blank sapph\\blank_sapph_mirror_VNATrc.002", 
                         "WR1o5 blank sapph\\blank_sapph_VNATrc.001", 
                         "WR1o5 blank sapph\\blank_sapph_VNATrc.002")

ITO_300_up = rt.reflm_std("WR1o5 ITO sapph\ITO_up_300_mirror_VNATrc.001", 
                          "WR1o5 ITO sapph\ITO_up_300_mirror_VNATrc.002", 
                          "WR1o5 ITO sapph\ITO_up_300_VNATrc.001", 
                          "WR1o5 ITO sapph\ITO_up_300_VNATrc.002")

ITO_300_down = rt.reflm_std("WR1o5 ITO sapph\ITO_down_300_mirror_VNATrc.001", 
                            "WR1o5 ITO sapph\ITO_down_300_mirror_VNATrc.002", 
                            "WR1o5 ITO sapph\ITO_down_300_VNATrc.001", 
                            "WR1o5 ITO sapph\ITO_down_300_VNATrc.002")

samps = {}
freqbd = {}
freqbd['wr3p4'] = [228e9,311e9]
freqbd['wr1p5'] = [490e9,680e9]

samps['ITO 150 nm up'] = rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.482e-3,window=100)
samps['ITO 150 nm up'].add_trc('wr3p4', ITO_150_up)


samps['ITO 150 nm down']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.482e-3,window=100)
samps['ITO 150 nm down'].add_trc('wr3p4', ITO_150_down)


samps['Blank, WR3.4']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.486e-3,window=100)
samps['Blank, WR3.4'].add_trc('wr3p4', blank_3o4)


samps['Blank, WR1.5']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.486e-3,window=100)
samps['Blank, WR1.5'].add_trc('wr1p5', blank_1o5)

samps['ITO 300 nm up']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.486e-3,window=100)
samps['ITO 300 nm up'].add_trc('wr1p5', ITO_300_up)

samps['ITO 300 nm down']=rt.trc(freqbounds=freqbd,n1=1,n2=1,d=.486e-3,window=100)
samps['ITO 300 nm down'].add_trc('wr1p5', ITO_300_down)


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


