import numpy as np
import matplotlib.pyplot as plt
import refltool2 as rt
import ncalclib as ncl
from os import chdir
from matplotlib.ticker import AutoMinorLocator

samps={}
freqbd={}
freqbd['wr1p5']=[490e9,680e9]
freqbd['wr2p2']=[330e9,490e9]
freqbd['wr3p4']=[228e9,311e9]
freqbd['wr5p1']=[170e9,228e9]
freqbd['wr6p5']=[110e9,170e9]
freqbd['wr10']=[70e9,110e9]


def f2c(mat):
    return (mat[:,0]+1j*mat[:,1]).reshape((len(mat[:,0]),1))

def c2f(mat):
    return np.hstack((mat.real, mat.imag))

def reflm_std(mirror_reference_filename,mirror_refl_filename,sample_reference_filename,sample_refl_filename):
    samp_normed=f2c(np.loadtxt(sample_refl_filename,skiprows=1)[:,3:5])/f2c(np.loadtxt(sample_reference_filename,skiprows=1)[:,3:5])
    mirr_normed=f2c(np.loadtxt(mirror_refl_filename,skiprows=1)[:,3:5])/f2c(np.loadtxt(mirror_reference_filename,skiprows=1)[:,3:5])
    reflp=c2f(samp_normed/mirr_normed);
    reflp[:,0]=-reflp[:,0]
    return np.hstack((
        np.loadtxt(sample_refl_filename,skiprows=1)[:,0:1], 
        0*np.loadtxt(sample_refl_filename,skiprows=1)[:,0:1],
        reflp))


bands = ["wr10", "wr6p5", "wr5p1", "wr3p4", "wr2p2", "wr1p5"]
for elem in bands:
    samps['blank' + elem]=rt.trc(freqbounds=freqbd, n1=1, n2=1, d=1, window=100)
    samps['gaas' + elem]=rt.trc(freqbounds=freqbd, n1=1, n2=1, d=1, window=100)
    samps['si_yellow' + elem]=rt.trc(freqbounds=freqbd, n1=1, n2=1, d=1, window=100)
    samps['si_pink' + elem]=rt.trc(freqbounds=freqbd, n1=1, n2=1, d=1, window=100)

for elem in bands:
    chdir("../" + elem)
    bnd = elem
    samps['blank' + elem].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "blank_VNATrc.001", "blank_VNATrc.002" ))
    samps['gaas' + elem].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "gaas_VNATrc.001", "gaas_VNATrc.002" ))
    samps['si_yellow' + elem].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_yellow_VNATrc.001", "si_yellow_VNATrc.002" ))
    samps['si_pink' + elem].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_pink_VNATrc.001", "si_pink_VNATrc.002" ))



'''
chdir("../wr10")
bnd="wr10"

samps['blank'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "blank_VNATrc.001", "blank_VNATrc.002" ))
samps['gaas'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "gaas_VNATrc.001", "gaas_VNATrc.002" ))
samps['si_yellow'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_yellow_VNATrc.001", "si_yellow_VNATrc.002" ))
samps['si_pink'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_pink_VNATrc.001", "si_pink_VNATrc.002" ))

chdir("../wr6p5")
bnd="wr6p5"

samps['blank'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "blank_VNATrc.001", "blank_VNATrc.002" ))
samps['gaas'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "gaas_VNATrc.001", "gaas_VNATrc.002" ))
samps['si_yellow'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_yellow_VNATrc.001", "si_yellow_VNATrc.002" ))
samps['si_pink'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_pink_VNATrc.001", "si_pink_VNATrc.002" ))

chdir("../wr5p1")
bnd="wr5p1"

samps['blank'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "blank_VNATrc.001", "blank_VNATrc.002" ))
samps['gaas'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "gaas_VNATrc.001", "gaas_VNATrc.002" ))
samps['si_yellow'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_yellow_VNATrc.001", "si_yellow_VNATrc.002" ))
samps['si_pink'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_pink_VNATrc.001", "si_pink_VNATrc.002" ))

chdir("../wr3p4")
bnd="wr3p4"

samps['blank'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "blank_VNATrc.001", "blank_VNATrc.002" ))
samps['gaas'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "gaas_VNATrc.001", "gaas_VNATrc.002" ))
samps['si_yellow'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_yellow_VNATrc.001", "si_yellow_VNATrc.002" ))
samps['si_pink'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_pink_VNATrc.001", "si_pink_VNATrc.002" ))

chdir("../wr2p2")
bnd="wr2p2"

samps['blank'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "blank_VNATrc.001", "blank_VNATrc.002" ))
samps['gaas'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "gaas_VNATrc.001", "gaas_VNATrc.002" ))
samps['si_yellow'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_yellow_VNATrc.001", "si_yellow_VNATrc.002" ))
samps['si_pink'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_pink_VNATrc.001", "si_pink_VNATrc.002" ))

chdir("../wr1p5")
bnd="wr1p5"

samps['blank'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "blank_VNATrc.001", "blank_VNATrc.002" ))
samps['gaas'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "gaas_VNATrc.001", "gaas_VNATrc.002" ))
samps['si_yellow'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_yellow_VNATrc.001", "si_yellow_VNATrc.002" ))
samps['si_pink'].add_trc(bnd, reflm_std("mirror_VNATrc.001", "mirror_VNATrc.002", "si_pink_VNATrc.001", "si_pink_VNATrc.002" ))
'''
chdir("../process")

plt.close('all')
#plt.plot(samps['si_yellow'].trc_bd_fstack()[:,0],samps['si_yellow'].bd_reflectivity()[:,2])
#plt.plot(samps['si_pink'].trc_bd_fstack()[:,0],samps['si_pink'].bd_reflectivity()[:,2])
for key in samps.keys():
    elem = samps[key]
    a = rt.savitzky_golay(elem.efield_surf_abs(), 51, 2, deriv=0, rate=1)
    plt.plot(elem.trc_bd_fstack()[:,0]*1e-9, a, label=key)
    np.savetxt(key, np.vstack((elem.trc_bd_fstack()[:,0]*1e-9, a)).T)

#a = rt.savitzky_golay(samps['blankwr1p5'].efield_surf_abs(), 51, 2, deriv=0, rate=1)
#b = rt.savitzky_golay(samps['gaaswr1p5'].efield_surf_abs(), 51, 2, deriv=0, rate=1)
#plt.plot(samps['blankwr1p5'].trc_bd_fstack()[:,0]/1e9,a)
#plt.plot(samps['gaaswr1p5'].trc_bd_fstack()[:,0]/1e9,b)
plt.ylim(0,2)
plt.show()

'''
fig, ax = plt.subplots()
plt.minorticks_on()
plt.plot(pellet1_1[:,0]/1e9,pellet1_1[:,2])
plt.plot(pellet1_2[:,0]/1e9,pellet1_2[:,2])
plt.plot(pellet1_3[:,0]/1e9,pellet1_3[:,2])
plt.plot(pellet1_4[:,0]/1e9,pellet1_4[:,2])
plt.plot(pellet1_5[:,0]/1e9,pellet1_5[:,2])
plt.ylim(-.5,.5)
plt.xlim(475,675)
plt.title("Reflectance, pellet 1")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Reflectance (Re)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig("apr7_pel1_re.pdf",format="pdf")
plt.clf()

fig, ax = plt.subplots()
plt.minorticks_on()
plt.plot(pellet1_1[:,0]/1e9,pellet1_1[:,3])
plt.plot(pellet1_2[:,0]/1e9,pellet1_2[:,3])
plt.plot(pellet1_3[:,0]/1e9,pellet1_3[:,3])
plt.plot(pellet1_4[:,0]/1e9,pellet1_4[:,3])
plt.plot(pellet1_5[:,0]/1e9,pellet1_5[:,3])
plt.ylim(-.5,.5)
plt.xlim(475,675)
plt.title("Reflectance, pellet 1")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Reflectance (Im)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig("apr7_pel1_im.pdf",format="pdf")
plt.clf()

fig, ax = plt.subplots()
plt.minorticks_on()
plt.plot(pellet2_1[:,0]/1e9,pellet2_1[:,2])
plt.plot(pellet2_2[:,0]/1e9,pellet2_2[:,2])
plt.plot(pellet2_3[:,0]/1e9,pellet2_3[:,2])
plt.plot(pellet2_4[:,0]/1e9,pellet2_4[:,2])
plt.plot(pellet2_5[:,0]/1e9,pellet2_5[:,2])
plt.ylim(-.5,.5)
plt.xlim(475,675)
plt.title("Reflectance, pellet 2")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Reflectance (Re)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig("apr7_pel2_re.pdf",format="pdf")
plt.clf()

fig, ax = plt.subplots()
plt.minorticks_on()
plt.plot(pellet2_1[:,0]/1e9,pellet2_1[:,3])
plt.plot(pellet2_2[:,0]/1e9,pellet2_2[:,3])
plt.plot(pellet2_3[:,0]/1e9,pellet2_3[:,3])
plt.plot(pellet2_4[:,0]/1e9,pellet2_4[:,3])
plt.plot(pellet2_5[:,0]/1e9,pellet2_5[:,3])
plt.ylim(-.5,.5)
plt.xlim(475,675)
plt.title("Reflectance, pellet 2")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Reflectance (Im)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig("apr7_pel2_im.pdf",format="pdf")
plt.clf()

fig, ax = plt.subplots()
plt.minorticks_on()
plt.plot(pellet3_1[:,0]/1e9,pellet3_1[:,2])
plt.plot(pellet3_2[:,0]/1e9,pellet3_2[:,2])
plt.plot(pellet3_3[:,0]/1e9,pellet3_3[:,2])
plt.plot(pellet3_4[:,0]/1e9,pellet3_4[:,2])
plt.plot(pellet3_5[:,0]/1e9,pellet3_5[:,2])
plt.ylim(-.5,.5)
plt.xlim(475,675)
plt.title("Reflectance, pellet 3")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Reflectance (Re)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig("apr7_pel3_re.pdf",format="pdf")
plt.clf()

fig, ax = plt.subplots()
plt.minorticks_on()
plt.plot(pellet3_1[:,0]/1e9,pellet3_1[:,3])
plt.plot(pellet3_2[:,0]/1e9,pellet3_2[:,3])
plt.plot(pellet3_3[:,0]/1e9,pellet3_3[:,3])
plt.plot(pellet3_4[:,0]/1e9,pellet3_4[:,3])
plt.plot(pellet3_5[:,0]/1e9,pellet3_5[:,3])
plt.ylim(-.5,.5)
plt.xlim(475,675)
plt.title("Reflectance, pellet 3")
plt.xlabel("Frequency (GHz)")
plt.ylabel("Reflectance (Im)")
ax.xaxis.set_minor_locator(AutoMinorLocator(2))
ax.yaxis.set_minor_locator(AutoMinorLocator(2))
plt.savefig("apr7_pel3_im.pdf",format="pdf")
plt.clf()
'''


