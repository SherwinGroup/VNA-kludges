import numpy as np

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
