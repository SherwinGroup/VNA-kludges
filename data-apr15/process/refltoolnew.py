import numpy as np
from scipy.optimize import minimize

'''
Things that need to be done here:
 - Make naming conventions for r and R consistent.  Complex reflection coefficient?
   Reflectivity?  Complex reflectivity?  Which one is power vs field?
 - At the beginning, make a list of variable names to clarify the above point
'''


########################################
#general utilities
########################################



def f2c(vector):
    """
    this will turn a numpy array vector of form [real, imag] to [complex]
    """
    return (vector[:, 0] + 1j * vector[:, 1]).reshape((len(vector[:, 0]), 1))

def c2f(vector):
    """
    The inverse of f2c
    """
    return np.hstack((vector.real, vector.imag))

#bwwuaaahhh
def load_refl_data(mirror_refer_fname, mirror_refl_fname, sample_refer_fname, 
              sample_refl_fname):
    """
    This function will load standard reflection data that has been normalized
    to a mirror reference.
    
    mirror_refer_fname: file name of the mirror reference data    
    mirror_refl_fname: file name of the mirror reflection data
    sample_refer_fname: file name of the sample reference data
    sample_refl_fname: file name of the mirror reflection data
    
    returns: I...don't...know...
    """
    sample_normed = f2c(np.loadtxt(sample_refl_fname)[:,3:5]) / f2c(np.loadtxt(sample_refer_fname)[:,3:5])
    mirror_normed = f2c(np.loadtxt(mirror_refl_fname)[:,3:5]) / f2c(np.loadtxt(mirror_refer_fname)[:,3:5])
    reflp = c2f(sample_normed / mirror_normed);
    reflp[:,0] = -reflp[:,0]
    return np.hstack((np.loadtxt(sample_refl_fname)[:,0:1], 
                      0 * np.loadtxt(sample_refl_fname)[:,0:1], reflp))

########################################
########################################


########################################
#refractive index calc. functions
########################################

def slab_refl_n(n, n1, n2, d, freq_range):
    """
    used to be called "slabr"
    
    calculates frequency-dependent complex reflection coefficient for a "slab" geometry
                  ->         |||                     |||
             n1   ->         |||          n          |||          n2
                  ->         |||                     |||
    d = thickness of middle layer
    freq range is desired frequency or array of frequencies
    
    returns: complex reflection coefficient back into material n1
    """
    r = ((n1-n)/(n1+n) + (n-n2)/(n+n2) * np.exp(4*np.pi*1j*n*d/2.9979e8*freq_range)) / 
        (1 + (n1-n)/(n1+n) * (n-n2)/(n+n2) * np.exp(4*np.pi*1j*n*d/2.9979e8*freq_range))
    return r

def slab_refl_r(n,r12,r23,d,freq):
    """
    used to be called "slabrr"
    
    This is the same as slab_refl_n, but using the "simple" complex reflection 
    coefficients
    """
    r = ((r12 + r23*np.exp(4*np.pi*1j*d*n*freq/2.9979e8)) / 
        (1 + r12*r23*np.exp(4*np.pi*1j*d*n*freq/2.9979e8)))
    return r
    
'''
calculate refractive index from reflectance
all units mks
'''

def index_errfunc(r_data, n, n1, n2, d, freq):
    """
    This takes complex reflection coefficient data and compares it to a calculated
    complex reflection coefficient from given n, n1 and n2.  I think this is to 
    provide a least-squares function for fitting to data?
    
    r_data: complex reflection coefficient from measurement
    (n, n1, n2, d, freq): imputs to the slab_refl_n function
    
    returns: the absolute square of the difference
    """
    diff_square = abs(slab_refl_n(n, n1, n2, d, freq) - r_data)**2
    return diff_square

def index_calc(r_data, n_guess, n1, n2, d, freq_range):
    """
    used to be called "ncalc".  This needs a lot of work, it's just about the
    least pythonic thing I've seen here, which is saying a lot.  Brayden 
    disapproves.
    
    This function will perform a fit to a set of reflection coefficient data
    that will approximate the complex index of refraction.  You can't invert 
    the slab_refl_n function, so you have to fit to it.  
    
    r_data: complex reflection coefficient from measurement
    n_guess: I have no idea what this is
    (n1, n2, d, freq_range): for slab_refl_n function
    
    returns: complex index of refraction
    """
    guess_array = np.array([n_guess.real, n_guess.imag])
    nout = np.array([])
    ii = 0
    if len(n1) != 1:
        print n1[ii]
        print guess_array
        for freq in freq_range:
            new_errfun = lambda x: index_errfunc(r_data[ii], x[0] + 1j*x[1], # Couldn't you just do f2c?
                                                 n1[ii, 0], n2, d, freq_range[ii])
            outp = minimize(new_errfun, guess_array, method='Nelder-Mead')
            nout = np.append(nout,[outp.x[0]+1j*outp.x[1]])
            ii += 1
    '''
    else:
        for freq in freq_range:
            new_errfun = lambda x: refr_errfun(r_data[ii],x[0]+1j*x[1],n1,n2,d,freq_range[ii])
            outp=minimize(new_errfun,guess_array,method='Nelder-Mead')
            nout=np.append(nout,[outp.x[0]+1j*outp.x[1]])
            ii=ii+1
    '''
    return nout

def sheet_refl(n1, n2, sigma):
    """
    This calculates the complex(?) reflection coefficient?
    
    n1: index of left material
    n2: index of right material
    sigma: surface conductivity at the interface
    Waves are travelling left-right
    
    returns: some sort of reflection coefficient
    """
    return (n2 - n1 + sigma/(2.9979e8*8.85e-12)) / (n2 + n1 + sigma/(2.9979e8*8.85e-12))

def sheet_trans(n1,n2,sigma):
    """
    This function calculates the complex transmittance
    
    n1: index of left material
    n2: index of right material
    sigma: surface conductivity at the interface
    Waves are travelling left-right
    
    returns: some sort of transmission coefficient
    """
    return np.sqrt(n2 / n1) * 2 * n1 / (n1 + n2 + sigma/(2.9979e8*8.85e-12))

#air on both sides
#slab of index n in between with surface sheets of conductance sigma
'''
What the fuck is the above comment for?
'''

def hr_rv(freq, sigma_dc):
    """
    No fucking clue
    """
    return 1 - 2 * np.sqrt(4 * 8.85e-12 * np.pi * freq / sigma_dc)

def sheet_slab_refl_n(n, sigma, d, freq):
    """
    Super confused here, too
    """
    r12 = sheet_trans(1, 1, sigma)
    return r12 * (1 - np.exp(4 * np.pi * 1j * n * d * freq/ 2.9979e8)) / 
           (1 - r12**2 * np.exp(4 * np.pi * 1j * n * d * freq / 2.9979e8))

def ncalc_back(r_data, n1, n_front, d_front, n2, d, freq_range, nguess):
    """
    This looks really confusing
    """
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

class trc(object): #means trace
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
def savitzky_golay(y, window_size, order, deriv=0, rate=1):
    """Smooth (and optionally differentiate) data with a Savitzky-Golay filter.
    The Savitzky-Golay filter removes high frequency noise from data.
    It has the advantage of preserving the original shape and
    features of the signal better than other types of filtering
    approaches, such as moving averages techniques.
    Parameters
    ----------
    y : array_like, shape (N,)
        the values of the time history of the signal.
    window_size : int
        the length of the window. Must be an odd integer number.
    order : int
        the order of the polynomial used in the filtering.
        Must be less then `window_size` - 1.
    deriv: int
        the order of the derivative to compute (default = 0 means only smoothing)
    Returns
    -------
    ys : ndarray, shape (N)
        the smoothed signal (or it's n-th derivative).
    Notes
    -----
    The Savitzky-Golay is a type of low-pass filter, particularly
    suited for smoothing noisy data. The main idea behind this
    approach is to make for each point a least-square fit with a
    polynomial of high order over a odd-sized window centered at
    the point.
    Examples
    --------
    t = np.linspace(-4, 4, 500)
    y = np.exp( -t**2 ) + np.random.normal(0, 0.05, t.shape)
    ysg = savitzky_golay(y, window_size=31, order=4)
    import matplotlib.pyplot as plt
    plt.plot(t, y, label='Noisy signal')
    plt.plot(t, np.exp(-t**2), 'k', lw=1.5, label='Original signal')
    plt.plot(t, ysg, 'r', label='Filtered signal')
    plt.legend()
    plt.show()
    References
    ----------
    .. [1] A. Savitzky, M. J. E. Golay, Smoothing and Differentiation of
       Data by Simplified Least Squares Procedures. Analytical
       Chemistry, 1964, 36 (8), pp 1627-1639.
    .. [2] Numerical Recipes 3rd Edition: The Art of Scientific Computing
       W.H. Press, S.A. Teukolsky, W.T. Vetterling, B.P. Flannery
       Cambridge University Press ISBN-13: 9780521880688
    """
    import numpy as np
    from math import factorial

    try:
        window_size = np.abs(np.int(window_size))
        order = np.abs(np.int(order))
    except ValueError, msg:
        raise ValueError("window_size and order have to be of type int")
    if window_size % 2 != 1 or window_size < 1:
        raise TypeError("window_size size must be a positive odd number")
    if window_size < order + 2:
        raise TypeError("window_size is too small for the polynomials order")
    order_range = range(order+1)
    half_window = (window_size -1) // 2
    # precompute coefficients
    b = np.mat([[k**i for i in order_range] for k in range(-half_window, half_window+1)])
    m = np.linalg.pinv(b).A[deriv] * rate**deriv * factorial(deriv)
    # pad the signal at the extremes with
    # values taken from the signal itself
    firstvals = y[0] - np.abs( y[1:half_window+1][::-1] - y[0] )
    lastvals = y[-1] + np.abs(y[-half_window-1:-1][::-1] - y[-1])
    y = np.concatenate((firstvals, y, lastvals))
    return np.convolve( m[::-1], y, mode='valid')