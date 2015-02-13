import numpy as np
from scipy.optimize import minimize

########################################
#general utilities
########################################

def f2c(mat):
    return (mat[:,0]+1j*mat[:,1]).reshape((len(mat[:,0]),1))

def c2f(mat):
    return np.hstack((mat.real, mat.imag))

#bwwuaaahhh
def reflm_std(mirror_reference_filename,mirror_refl_filename,sample_reference_filename,sample_refl_filename):
    samp_normed=f2c(np.loadtxt(sample_refl_filename)[:,3:5])/f2c(np.loadtxt(sample_reference_filename)[:,3:5])
    mirr_normed=f2c(np.loadtxt(mirror_refl_filename)[:,3:5])/f2c(np.loadtxt(mirror_reference_filename)[:,3:5])
    reflp=c2f(samp_normed/mirr_normed);
    reflp[:,0]=-reflp[:,0]
    return np.hstack((
        np.loadtxt(sample_refl_filename)[:,0:1], 
        0*np.loadtxt(sample_refl_filename)[:,0:1],
        reflp))

########################################
########################################


########################################
#refractive index calc. functions
########################################

def slabr(n,n1,n2,d,freq_range):
    r=((n1-n)/(n1+n)+(n-n2)/(n+n2)*np.exp(4*np.pi*1j*n*d/2.9979e8*freq_range))/(1+(n1-n)/(n1+n)*(n-n2)/(n+n2)*np.exp(4*np.pi*1j*n*d/2.9979e8*freq_range))
    return r

def slabrr(n,r12,r23,d,freq):
    r=((r12+r23*np.exp(4*np.pi*1j*d*n*freq/2.9979e8))/(1+r12*r23*np.exp(4*np.pi*1j*d*n*freq/2.9979e8)))

'''
calculate refractive index from reflectance
all units mks
'''

def refr_errfun(r_data,n,n1,n2,d,freq):
    out=abs(slabr(n,n1,n2,d,freq)-r_data)**2
    return out

def ncalc(r_data,nguess,n1,n2,d,freq_range):
    guess_array=np.array([nguess.real, nguess.imag])
    nout=np.array([])
    ii=0
    if len(n1) != 1:
        print n1[ii]
        print guess_array
        for freq in freq_range:
            new_errfun = lambda x: refr_errfun(r_data[ii],x[0]+1j*x[1],n1[ii,0],n2,d,freq_range[ii])
            outp=minimize(new_errfun,guess_array,method='Nelder-Mead')
            nout=np.append(nout,[outp.x[0]+1j*outp.x[1]])
            ii=ii+1
    '''
    else:
        for freq in freq_range:
            new_errfun = lambda x: refr_errfun(r_data[ii],x[0]+1j*x[1],n1,n2,d,freq_range[ii])
            outp=minimize(new_errfun,guess_array,method='Nelder-Mead')
            nout=np.append(nout,[outp.x[0]+1j*outp.x[1]])
            ii=ii+1
    '''
    return nout

def sheet_refl(n1,n2,sigma):
    #sigma is the surface conductivity
    #freqs is a vector
    return (n2+sigma/2.9979e8/8.85e-12-n1)/(n2+sigma/2.9979e8/8.85e-12+n1)

def sheet_trans(n1,n2,sigma):
    #transmittance
    return np.sqrt(n2/n1)*2*n1/(n1+n2+sigma/2.9979e8/8.85e-12)

#air on both sides
#slab of index n in between with surface sheets of conductance sigma

def hr_rv(freq,sigma_dc):
    return (1- 2*np.sqrt(4*8.85e-12*np.pi*freq/sigma_dc))

def slabr_sheet(n,sigma,d,freq):
    r12=sheet_trans(1,1,sigma)
    return (r12*(1-np.exp(4*np.pi*1j*n*d/2.9979e8*freq))/
            (1-r12**2*np.exp(4*np.pi*1j*n*d/2.9979e8*freq)))

def ncalc_back(r_data,n1,n_front,d_front,n2,d,freq_range,nguess):
    r12_front=(n1-n_front)/(n1+n_front)
    guess_array=np.array([nguess.real, nguess.imag])
    nout=np.array([])
    ii=0
    for freq in freq_range:
        new_errfun = lambda x: (
            abs(
                (r12_front[ii]+slabr(x[0]+x[1]*1j,n_front[ii],n2,d,freq)*np.exp(
                    4*np.pi*1j*n_front[ii]*d_front*freq/2.9979e8
                    )
                )/
                (1+r12_front[ii]*slabr(x[0]+x[1]*1j,n_front[ii],n2,d,freq)*np.exp(
                    4*np.pi*1j*n_front[ii]*d_front*freq/2.9979e8
                    )
                 ) - r_data[ii])**2
        )
        outp=minimize(new_errfun,guess_array,method='Nelder-Mead')
        nout=np.append(nout,[outp.x[0]+1j*outp.x[1]])
        ii=ii+1
    return nout

########################################
########################################


########################################
#drude-related functions
########################################

#calculate a drude dielectric function for a set of parameters
#inputs are freqs -- the frequencies at which the dielectric function is calculated
#freqs should be an array
#fp is the plasma frequency
#gamma is 1/tau -- the relaxation rate

def drude_eps(freqs,fp,gamma):
    return ((1-4*np.pi**2*fp**2/(4*np.pi**2*freqs**2+gamma**2))+1j*(4*np.pi**2*gamma*fp**2/(8*np.pi**3*freqs**3+4*np.pi**2*freqs**2*gamma))).reshape((len(freqs),1))

def drude_eps(freqs,fp,gamma):
    return ((1-4*np.pi**2*fp**2/(4*np.pi**2*freqs**2+gamma**2))+1j*(4*np.pi**2*gamma*fp**2/(8*np.pi**3*freqs**3+4*np.pi**2*freqs**2*gamma))).reshape((len(freqs),1)).view(float)


#freqs is a 1d array of freqs
# eps_in is a complex 1d array of measured dielectric function
def drude_eps_errfun(freqs,eps_in,fp,gamma):
    return np.trapz(np.abs(drude_eps(freqs,fp,gamma)-eps_in)**2,axis=0)

#spec should be an array in the same form as the above
def drude_spec_fit(freqs,eps_in,fp_guess,gamma_guess):
    from scipy.optimize import minimize
    #input_array is of the form: [wp_guess gamma_guess]
    def testf(input_array):
        return drude_eps_errfun(freqs,eps_in,input_array[0],input_array[1])
    return minimize(testf,np.array([fp_guess, gamma_guess]),method='Nelder-Mead',tol=1)

def drude_spec_fit2(freqs,eps_in,fp_guess,gamma_guess):
    from scipy.optimize import curve_fit
    return curve_fit(drude_eps2,freqs,eps_in)

########################################
########################################


########################################
#main part:
#class for reflectance data
########################################

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
        try:
            return self.trc[band][self.bd_ind(band)[0]:self.bd_ind(band)[1]]
        except ValueError:
            return self.trc[band]
		
    def list_bands(self):
        outp=[]
        for band in self.trc:
            outp.append(band)
        return outp 
		
    def trc_bd_fstack(self):
        minfreqlist={}
        for band in self.list_bands():
            minfreqlist[band]=self.bdfreq(band)[0]
        try:
            outp=self.bdtrc(min(minfreqlist,key=minfreqlist.get))
            minfreqlist.pop(min(minfreqlist,key=minfreqlist.get),None)
            while minfreqlist != {}:
                nextband=min(minfreqlist,key=minfreqlist.get) 
                outp=np.vstack((outp,self.bdtrc(nextband)))
                minfreqlist.pop(nextband,None)
            return outp[outp[:,0].argsort()]
        except ValueError:
            bandlist=self.list_bands()
            outp=self.trc[bandlist[0]]
            bandlist.pop(0)
            while len(bandlist) != 0:
                outp=np.vstack((outp,self.trc[bandlist(0)]))
                bandlist.pop(0)
            return outp[out[:,0].argsort()]

    def nguess_make(self):
        from scipy.signal import argrelmax
        peakind=argrelmax(np.abs(self.trc_bd_fstack().view(complex)[:,1])**2,order=self.window)
        #print peakind
        return (2.9979e8/np.diff(self.trc_bd_fstack()[peakind,0]).mean(axis=1)/2/self.d)[0]
        
    def n_slab(self):
        return ncalc(self.trc_bd_fstack().view(complex)[:,1],complex(self.nguess_make()),self.n1,self.n2,self.d,self.trc_bd_fstack()[:,0]).view(float).reshape((len(self.trc_bd_fstack()[:,0]),2))
		
    def n_slab_manual_guess(self,nguess,kguess):
        return ncalc(self.trc_bd_fstack().view(complex)[:,1],nguess+1j*kguess,self.n1,self.n2,self.d,self.trc_bd_fstack()[:,0]).view(float).reshape((len(self.trc_bd_fstack()[:,0]),2))     
 
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

    def efield_surf_abs(self):
        return abs(1+f2c(self.trc_bd_fstack()[:,2:4])[:,0])

    def bd_reflectivity(self):
        return np.hstack((
            self.trc_bd_fstack()[:,0:2],
            np.abs(
                self.trc_bd_fstack()[:,4:5]**2*4*self.trc_bd_fstack()[:,2:3]**2+
                self.trc_bd_fstack()[:,5:6]**2*4*self.trc_bd_fstack()[:,3:4]**2)
            ))

    def hr_fit_raw(self):
        from scipy.optimize import curve_fit
        return curve_fit(hr_rv,self.trc_bd_fstack()[:,0],self.bd_reflectivity()[:,2])
        
    def bd_rv_slab2inf(self):
        return np.abs((1-self.n_slab().view(complex))/(1+self.n_slab().view(complex)))
        
#for more complicated meas: add fmt tuple for handling multiple measurements of mirror, sample or reference
#for now only averages reflection coefficients
#filename should be formatted as a float
#def trc_fromblk(filename,**kwargs):

########################################
########################################
