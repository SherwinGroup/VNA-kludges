import numpy as np
import refltool as rt
from pylab import *

freq=np.linspace(0e9,1e12,num=1000)

n_ito_1=5+5j
n_ito_2=50+50j

n1=3.
n2=1.

d1=150e-9
d2=300e-9
d3=450e-9

outp1=rt.slabr(n_ito_2,n1,n2,d1,freq)
outp2=rt.slabr(n_ito_2,n1,n2,d2,freq)

outp3=rt.slabr(n_ito_2,n1,n2,d3,freq)

'''
plot(freq,outp1.real)
plot(freq,outp1.imag)
plot(freq,outp2.real)
plot(freq,outp2.imag)
plot(freq,outp3.real)
plot(freq,outp3.imag)
show()
'''

freqs=500e9
d=np.linspace(100e-9,1000e-9,num=1000)

dsaph1=1e-3*3.7
dsaph2=1e-3*3.75
dsaph3=1e-3*3.98

outp=rt.slabr(n_ito_2,n1,n2,d,freqs)

#neff=(n2-outp)/(n2+outp)

r23=outp

'''
plot(d*1e9,outp.real)
plot(d*1e9,outp.imag)
show()
'''

r12=(1.-n1)/(1.+n1)

beta1=2*np.pi*n1*dsaph1*freqs/3e8
beta2=2*np.pi*n1*dsaph2*freqs/3e8
beta3=2*np.pi*n1*dsaph3*freqs/3e8

refl1=(r12+r23*np.exp(2.*1j*beta1))/(1+r12*r23*np.exp(2.*1j*beta1))
refl2=(r12+r23*np.exp(2.*1j*beta2))/(1+r12*r23*np.exp(2.*1j*beta2))
refl3=(r12+r23*np.exp(2.*1j*beta3))/(1+r12*r23*np.exp(2.*1j*beta3))

#plot(d*1e9,refl.real)
#plot(d*1e9,refl.imag)

#plot(d*1e9,abs(refl1+1))
plot(d*1e9,abs(refl2+1))
#plot(d*1e9,abs(refl3+1))
show()
