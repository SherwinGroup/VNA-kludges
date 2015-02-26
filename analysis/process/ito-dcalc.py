import numpy as np
import refltool as rt
import matplotlib.pyplot as plt



n1 = 3.
n2 = 1.
n_ito_1 = 67+67j
n_ito_2 = 100+100j

freq_array = np.linspace(1e9, 1e12, num=1000)
dist_array = np.linspace(100e-9, 500e-9, num=1000)
'''
outp1=rt.slabr(n_ito_2,n1,n2,d1,freq)
outp2=rt.slabr(n_ito_2,n1,n2,d2,freq)

outp3=rt.slabr(n_ito_2,n1,n2,d3,freq)

plot(freq,outp1.real)
plot(freq,outp1.imag)
plot(freq,outp2.real)
plot(freq,outp2.imag)
plot(freq,outp3.real)
plot(freq,outp3.imag)
show()
'''
single_dist = 150e-9
single_freq = 510e9

sapph_thickness = 495e-6

outp1=rt.slabr(n_ito_1, n1, n2, single_dist, freq_array)
outp2=rt.slabr(n_ito_2 ,n1, n2, single_dist, freq_array)
#neff=(n2-outp)/(n2+outp)

r23_1 = outp1
r23_2 = outp2
'''
plot(d*1e9,outp.real)
plot(d*1e9,outp.imag)
show()
'''

r12=(1.-n1)/(1.+n1)

#beta1=2*np.pi*n1*dsaph1*freqs/3e8
beta2=2*np.pi*n1*sapph_thickness*freq_array/3e8
#beta3=2*np.pi*n1*dsaph3*freqs/3e8

refl1=(r12+r23_1*np.exp(2.*1j*beta2))/(1+r12*r23_1*np.exp(2.*1j*beta2))
refl2=(r12+r23_2*np.exp(2.*1j*beta2))/(1+r12*r23_2*np.exp(2.*1j*beta2))

#plot(d*1e9,refl.real)
#plot(d*1e9,refl.imag)
plt.close()
#plt.plot(dist_array*1e9,abs(refl1+1))
plt.plot(freq_array*1e-9,abs(refl1+1))
#plot(d*1e9,abs(refl3+1))
plt.show()
