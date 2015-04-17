import numpy as np
import ncalclib as ncl
import drudefit as drude


class trc(object):
    def __init__(self,**kwargs):
        self.trc={}
        if ("freqbounds" in kwargs):
            self.freqbounds=kwargs['freqbounds']
        if ("n1" in kwargs):
            self.n1=kwargs['n1']
        if ("n2" in kwargs):
            self.n2=kwargs['n2']
        if ("d" in kwargs):
            self.d=kwargs['d']
        if ("window" in kwargs):
            self.window=kwargs['window']
    
    #format for data input into add_trc is
    #(Re(freq)) (Im(freq)) (Re(refl)) (Im(refl)) (std(Re(refl))) (std(Im(refl)))

    def add_trc(self,band,data):
        self.trc[band]=data
		
    def freq(self,band):
        return self.trc[band][:,0]
		
    def refl(self,band):
        return self.trc[band][:,2:4]
		
    def refl_std(self,band):
        return self.trc[band][:,4:6]
		
    def bd_ind(self,band):
        try:
            return [np.where(self.freq(band)<=self.freqbounds[band][0])[0][-1], np.where(self.freq(band)>=self.freqbounds[band][1])[0][0]]
        except IndexError:
            return [0, len(self.freq(band))]
		
    def bdfreq(self,band):
            return self.freq(band)[self.bd_ind(band)[0]:self.bd_ind(band)[1]]

		
    def bdrefl(self,band):
        return self.refl(band)[self.bd_ind(band)[0]:self.bd_ind(band,self)[1]]
		
    #attempts to match together entire spectrum
    def bdtrc(self,band):
        return self.trc[band][self.bd_ind(band)[0]:self.bd_ind(band)[1]]
		
    def list_bands(self):
        outp=[]
        for band in self.trc:
            outp.append(band)
        return outp 
		
    def trc_bd_fstack(self):
        minfreqlist={}
        for band in self.list_bands():
            minfreqlist[band]=self.bdfreq(band)[0]
        if minfreqlist != {}:
            outp=self.bdtrc(min(minfreqlist,key=minfreqlist.get))
            minfreqlist.pop(min(minfreqlist,key=minfreqlist.get),None)
            while minfreqlist != {}:
                nextband=min(minfreqlist,key=minfreqlist.get) 
                outp=np.vstack((outp,self.bdtrc(nextband)))
                minfreqlist.pop(nextband,None)
        return outp[outp[:,0].argsort()]
		
    def nguess_make(self):
        from scipy.signal import argrelmax
        peakind=argrelmax(np.abs(self.trc_bd_fstack().view(complex)[:,1])**2,order=self.window)
        #print peakind
        return (2.9979e8/np.diff(self.trc_bd_fstack()[peakind,0]).mean(axis=1)/2/self.d)[0]
        
    def n_slab(self):
        return ncl.ncalc(self.trc_bd_fstack().view(complex)[:,1],complex(self.nguess_make()),self.n1,self.n2,self.d,self.trc_bd_fstack()[:,0]).view(float).reshape((len(self.trc_bd_fstack()[:,0]),2))
		
    def n_slab_manual_guess(self,nguess,kguess):
        return ncl.ncalc(self.trc_bd_fstack().view(complex)[:,1],nguess+1j*kguess,self.n1,self.n2,self.d,self.trc_bd_fstack()[:,0]).view(float).reshape((len(self.trc_bd_fstack()[:,0]),2))     
 
    def n_inf(self):
        return np.hstack((
        np.real((1-self.trc_bd_fstack().view(complex)[:,1])/(1+self.trc_bd_fstack().view(complex)[:,1])).reshape((len(self.trc_bd_fstack()[:,0]),1)), 
        np.imag((1-self.trc_bd_fstack().view(complex)[:,1])/(1+self.trc_bd_fstack().view(complex)[:,1])).reshape((len(self.trc_bd_fstack()[:,0]),1))
        ))

    def n_inf_std(self):
        rp=self.trc_bd_fstack()[:,2]
        rpp=self.trc_bd_fstack()[:,3]
        sp=self.trc_bd_fstack()[:,4]
        spp=self.trc_bd_fstack()[:,5]

        denom=(rpp**2+rp**2+2*rp+1)**2
        
        np_rp=2*(rpp-rp-1)*(rpp+rp+1)/denom
        np_rpp=-4*(rp+1)*rpp/denom
        npp_rp=4*(rp+1)*rpp/denom
        npp_rpp=2*(rpp-rp-1)*(rpp+rp+1)/denom
        return np.hstack((
            np.sqrt(np_rp**2*sp**2+np_rpp**2*spp**2).reshape((len(self.trc_bd_fstack()[:,0]),1)),
            np.sqrt(npp_rp**2*sp**2+npp_rpp**2*spp**2).reshape((len(self.trc_bd_fstack()[:,0]),1))
        ))
		
    def eps_slab(self):
        return (self.n_slab().view(complex)**2).view(float)
		
    def epsq_slab(self):
        return np.hstack((np.polyval(np.polyfit(self.trc_bd_fstack()[:,0],self.eps_slab()[:,0],2),self.trc_bd_fstack()[:,0]).reshape((len(self.trc_bd_fstack()[:,0]),1)), np.polyval(np.polyfit(self.trc_bd_fstack()[:,0],self.eps_slab()[:,1],2),self.trc_bd_fstack()[:,0]).reshape((len(self.trc_bd_fstack()[:,0]),1))))

		
    def eps_inf(self):
        return (self.n_inf().view(complex)**2).view(float)

    def sig_ohm_inf(self):
        return np.imag(self.eps_inf().view(complex)*8.85e-12*2*np.pi*self.trc_bd_fstack()[:,0:1]).reshape(((len(self.trc_bd_fstack()[:,0]),1)))


    def bd_reflectivity(self):
        return np.hstack((
            self.trc_bd_fstack()[:,0:2],
            np.abs(
                    self.trc_bd_fstack()[:,2]**2 + self.trc_bd_fstack()[:,3]**2
            ).reshape((
                len(self.trc_bd_fstack()[:,2]), 1
                ))
        ))

    def hr_fit_raw(self):
        from scipy.optimize import curve_fit
        return curve_fit(ncl.hr_rv,self.trc_bd_fstack()[:,0],self.bd_reflectivity()[:,2])

    def hr_rho_raw(self):
        from scipy.optimize import curve_fit
        return curve_fit(ncl.hr_rho, self.trc_bd_fstack()[:,0],self.bd_reflectivity()[:,2])

    def hr_rho_muohm_raw(self):
        from scipy.optimize import curve_fit
        return curve_fit(ncl.hr_rho_muohm, self.trc_bd_fstack()[:,0],self.bd_reflectivity()[:,2])
        
    def bd_rv_slab2inf(self):
        return np.abs((1-self.n_slab().view(complex))/(1+self.n_slab().view(complex)))
        

#for more complicated meas: add fmt tuple for handling multiple measurements of mirror, sample or reference
#for now only averages reflection coefficients
#filename should be formatted as a float
#def trc_fromblk(filename,**kwargs):
