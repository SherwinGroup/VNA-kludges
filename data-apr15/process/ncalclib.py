import numpy as np
from scipy.optimize import minimize


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

def hr_rho(freq,rho_dc): #rho in ohm-cm
    return (1-.2*np.sqrt(4*8.85e-12*np.pi*freq*rho_dc))

def hr_rho_muohm(freq,rho_dc): #rho in microohm-cm
    return (1-.0002*np.sqrt(4*8.85e-12*np.pi*freq*rho_dc))

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


#smart wrapper for ncalc
#calculates full n spectrum then alters guess to converge on a particular branch
'''
def nsm
'''
