import numpy as np
import matplotlib.pyplot as plt
import refltool as rt

'''
Calculate transfer matrix for an input list of optical parameters at a particular frequency.
parlist is an np array of the form

[[ 0                            eps of start medium    mu of start medium    ]
 [ thickness of first medium    eps of first medium    mu of first medium    ]
 [ thickness of second medium   eps of second medium   mu of second medium   ]
   ...
 [ thickness of nth medium      eps of nth medium      mu of nth medium      ]
 [ 0                            eps of end medium      mu of end medium      ]].
 
epsilon is the relative permittivity
mu is the relative permeability

Both mu and epsilon may be complex!
'''

'''
This one does not calculate it as a function of distance
'''
def TransferMatrix(ParList, freq):

	TransfMatrix = np.array([[1., 0.], [0., 1.]])
	
	for Slab_i in np.arange(ParList.shape[0]-1):
		#print Slab_i
		n_slab = np.sqrt(ParList[Slab_i, 1] * ParList[Slab_i, 2]) #refractive index inside material
		
		'''
		PropagatorMatrix corresponds to calculating the electric field at the far end of the slab
		knowing the electric field at the front of the slab.
		The form of the matrix is
		[[ e^(-i k l) 0        ]
		 [ 0          e^(i k l)]], 
		where l is the thickness of the slab and k is the wave vector inside the slab.
		k is related to the free-space wave vector k0 by k = k0 * n,
		where n is the refractive index inside the slab:
		n = sqrt(mu * epsilon).
		'''
		PropagatorMatrix = np.array([  [np.exp(-1j * 4.*np.pi * ParList[Slab_i, 0.] * n_slab * freq / 2.9979e8), 0.], 
			[0., np.exp(1j * 4.*np.pi * ParList[Slab_i, 0] * n_slab * freq / 2.9979e8)]  ])
		
		'''
		InterfaceMatrix uses the boundary conditions to calculate the fields in the next interface from the fields
		at the current interface.
		The form is
		1/t * [[ 1 r ]
		       [ r 1 ]].
		Here, t and r are the Fresnel transmission and reflection coefficients associated with the interface.
		r is given by (Z2 - Z1) / (Z2 + Z1)
		and t is given by  2*Z2 / (Z2 + Z1),
		where Z2 is the impedance of the next slab and Z1 is the impedance of the current slab.
		The impedance Z of a medium is given by Z = np.sqrt( mu / eps).
		'''
		Z1 = np.sqrt(ParList[Slab_i,     2]/ParList[Slab_i,     1])
		Z2 = np.sqrt(ParList[Slab_i + 1, 2]/ParList[Slab_i + 1, 1])
		
		t = 2.*Z2 /(Z2 + Z1)
		r = (Z2 - Z1)/(Z2 + Z1)
		
		InterfaceMatrix = 1/t * np.array([ [1., r], [r, 1.] ])
		
		'''
		Multiply the matrices:
		Propagator * Interface * Propagator * Interface * ....
		'''
		
		TransfMatrix = np.dot(TransfMatrix, np.dot(PropagatorMatrix, InterfaceMatrix))
	return TransfMatrix


'''
This one attempts to calculate the transfer matrix up to a particular distance TargetDistance along the stack of media
'''
def TransferMatrix_Before(ParList, TargetDistance, freq):
	#Initialize transfer matrix
	TransfMatrix_Before = np.array([ [1., 0.], [0., 1.] ])
	SlabBoundaries = ParList[:,0]
	
	if TargetDistance < 0:
		SlabFlag = 0
		TransfMatrix_Before = np.array([  [np.exp(-1j * 4.*np.pi * TargetDistance * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8), 0.], 
				[0., np.exp(1j * 4.*np.pi * TargetDistance * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8)]  ]),
	elif TargetDistance >= 0:
		'''
		The following loop determines which slab the target distance falls in.
		'''
		SlabFlag = -1
		TotalDistance = 0
		for Slab_Index in np.arange(len(SlabBoundaries)):
			if SlabFlag == -1:
				TotalDistance = TotalDistance + SlabBoundaries[Slab_Index]
				if TargetDistance - TotalDistance < 0:
					SlabFlag = Slab_Index
		'''
		if SlabFlag is still -1, then return the whole transfer matrix plus whatever extra distance there might be
		'''
		if SlabFlag == -1:
			TransfMatrix_Before = np.dot(
				TransferMatrix(ParList, freq),
				np.array([  [np.exp(-1j * 4.*np.pi * (TargetDistance - TotalDistance) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8), 0.], 
					[0., np.exp(1j * 4.*np.pi * (TargetDistance - TotalDistance) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8)]  ])
			)
		else: 
			ParList_Truncated = np.copy(ParList[0:(SlabFlag + 1),:])
			TransfMatrix_Before = np.dot(
				TransferMatrix(ParList_Truncated, freq),
				np.array([  [np.exp(-1j * 4.*np.pi * (TargetDistance - TotalDistance - SlabBoundaries[Slab_Index]) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8), 0.], 
					[0., np.exp(1j * 4.*np.pi * (TargetDistance - TotalDistance - SlabBoundaries[Slab_Index]) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8)]  ])
			)	
	return TransfMatrix_Before

	
'''
This one attempts to calculate the transfer matrix after a particular distance TargetDistance along the stack of media
'''
def TransferMatrix_After(ParList, TargetDistance, freq):
	#Initialize transfer matrix
	TransfMatrix_After = np.array([ [1., 0.], [0., 1.] ])
	SlabBoundaries = ParList[:,0]
	
	if TargetDistance < 0:
		SlabFlag = 0
		TransfMatrix_After = np.dot(
			np.array([  [np.exp(1j * 4.*np.pi * TargetDistance * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8), 0.], 
				[0., np.exp(-1j * 4.*np.pi * TargetDistance * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8)]  ]),
			TransferMatrix(ParList, freq)
		)
	elif TargetDistance >= 0:
		'''
		The following loop determines which slab the target distance falls in.
		'''
		SlabFlag = -1
		TotalDistance = 0
		for Slab_Index in np.arange(len(SlabBoundaries)):
			if SlabFlag == -1:
				TotalDistance = TotalDistance + SlabBoundaries[Slab_Index]
				if TargetDistance - TotalDistance < 0:
					SlabFlag = Slab_Index
		'''
		if SlabFlag is still -1, then return the whole transfer matrix plus whatever extra distance there might be.
		'''
		if SlabFlag == -1:
			TransfMatrix_After = np.array([  [np.exp(1j * 4.*np.pi * (TargetDistance - TotalDistance ) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8), 0.], 
				[0., np.exp(-1j * 4.*np.pi * (TargetDistance - TotalDistance) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8)]  ])
		else: 
			ParList_Truncated = np.copy(ParList[SlabFlag:len(SlabBoundaries),:])
			ParList_Truncated[0,0] = 0.
			TransfMatrix_After = np.dot(
				np.array([  [np.exp(1j * 4.*np.pi * (TargetDistance - TotalDistance) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8), 0.], 
					[0., np.exp(-1j * 4.*np.pi * (TargetDistance - TotalDistance) * np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * freq / 2.9979e8)]  ]),
				TransferMatrix(ParList_Truncated, freq)
			)
	
	return TransfMatrix_After

'''
To calculate the reflection coefficient from the transfer matrix, divide entry
(1, 0) by entry (0, 0).
'''

'''
EFieldCalc calculates the electric field along the length of the stack of materials.
ParList is the same array specifying the optical parameters of the stack.
freq is the frequency at which the electric field is to be calculated.
numpts is the number of points along the length of the stack at which the field is to be calculated.

EFieldCalc assumes that the incident electric field = 1; this function will calculate the factor by which
the incident electric field is multiplied.
'''
def EFieldCalc(ParList, TargetDistance, freq):
	TransfMatrix = TransferMatrix(ParList, freq)
	TransfMatrix_After = TransferMatrix_After(ParList, TargetDistance, freq)
	EFieldVec = np.dot(TransfMatrix_After, np.array([ [1./TransfMatrix[0, 0]], [0.] ]) )
	
	return EFieldVec[0, 0] + EFieldVec[1, 0]
	
def BFieldCalc(ParList, TargetDistance, freq):
	TransfMatrix = TransferMatrix(ParList, freq)
	TransfMatrix_After = TransferMatrix_After(ParList, TargetDistance, freq)
	
	SlabBoundaries = ParList[:,0]
	SlabFlag = -1
	TotalDistance = 0
	for Slab_Index in np.arange(len(SlabBoundaries)):
		if SlabFlag == -1:
			TotalDistance = TotalDistance + SlabBoundaries[Slab_Index]
			if TargetDistance - TotalDistance < 0:
				SlabFlag = Slab_Index
	
	
	return np.sqrt(ParList[SlabFlag, 1] * ParList[SlabFlag, 2]) * EFieldCalc(ParList, TargetDistance, freq)
	
	
'''
For the auxiliary fields, assume that the media are isotropic.
'''
def DFieldCalc(ParList, TargetDistance, freq):
	TransfMatrix = TransferMatrix(ParList, freq)
	TransfMatrix_After = TransferMatrix_After(ParList, TargetDistance, freq)
	EFieldVec = np.dot(TransfMatrix_After, np.array([ [1./TransfMatrix[0, 0]], [0.] ]) )
	
	SlabBoundaries = ParList[:,0]
	SlabFlag = -1
	TotalDistance = 0
	for Slab_Index in np.arange(len(SlabBoundaries)):
		if SlabFlag == -1:
			TotalDistance = TotalDistance + SlabBoundaries[Slab_Index]
			if TargetDistance - TotalDistance < 0:
				SlabFlag = Slab_Index
	
	return ParList[SlabFlag, 1] * (EFieldVec[0, 0] + EFieldVec[1, 0])
	
	
def HFieldCalc(ParList, TargetDistance, freq):
	TransfMatrix = TransferMatrix(ParList, freq)
	TransfMatrix_After = TransferMatrix_After(ParList, TargetDistance, freq)
	
	SlabBoundaries = ParList[:,0]
	SlabFlag = -1
	TotalDistance = 0
	for Slab_Index in np.arange(len(SlabBoundaries)):
		if SlabFlag == -1:
			TotalDistance = TotalDistance + SlabBoundaries[Slab_Index]
			if TargetDistance - TotalDistance < 0:
				SlabFlag = Slab_Index
	
	
	return np.sqrt(ParList[SlabFlag, 2] / ParList[SlabFlag, 1]) * EFieldCalc( ParList, TargetDistance, freq)




#example calculation


FreqList = np.linspace(235e9, 245e9, num=100)
ReflList = np.zeros(len(FreqList),dtype = np.complex128)

#d_critical = 3e8/240e9/4./3.4
d_critical = 1e-3

'''
ParList = np.array([
	[0., 1., 1.],
	[1e-3, 1., 1.],
	[d_critical, 3.4, 1.],
	[1e-4, 1., 1.], 
	[1e-4, (1+1j)/np.sqrt(2)*1e5, 1.],
	[1e-3, 1., 1.],
	[0., 1., 1.]
])
'''

ParList = np.array([
	[0., 1., 1.],
	[1e-3, 1., 1.],
	[1e-3, 11.7, 1.],
	[.5e-3, 11.7, 14],
	[1e-3, 1., 1.],
	[0., 1., 1.]
])

'''
#Freq = 240e9
#print TransferMatrix(ParList, Freq)



for freq_i in np.arange(len(FreqList)):
	Transf_i = TransferMatrix(ParList, FreqList[freq_i])
	ReflList[freq_i] = Transf_i[1, 0]/Transf_i[0, 0]

print Transf_i.dtype

plt.plot(FreqList, ReflList.real)
plt.plot(FreqList, ReflList.imag)
plt.show()
'''

#print TransferMatrix_Before(ParList, .1e-3, 240e9)
#print TransferMatrix_After(ParList, .5e-3, 240e9)

#print np.dot(TransferMatrix_Before(ParList, .5e-3, 240e9), TransferMatrix_After(ParList, .5e-3, 240e9))
#print TransferMatrix(ParList, 240e9)

Freq = 240e9
#DistList = np.linspace(1e-3+d_critical-.01e-3, 1e-3+d_critical+.01e-3, num=1000)
DistList = np.linspace(0,3e-3, num=1000)
EFieldList = np.zeros(len(DistList),dtype = np.complex128)
BFieldList = np.zeros(len(DistList),dtype = np.complex128)
DFieldList = np.zeros(len(DistList),dtype = np.complex128)
HFieldList = np.zeros(len(DistList),dtype = np.complex128)

for dist_i in np.arange(len(DistList)):
	EFieldList[dist_i] = EFieldCalc(ParList, DistList[dist_i], Freq)
	BFieldList[dist_i] = BFieldCalc(ParList, DistList[dist_i], Freq)
	DFieldList[dist_i] = DFieldCalc(ParList, DistList[dist_i], Freq)
	HFieldList[dist_i] = HFieldCalc(ParList, DistList[dist_i], Freq)


'''
print TransferMatrix_After(ParList, .000000001, Freq)
print TransferMatrix(ParList, Freq)

print TransferMatrix_After(ParList, 1.5e-3, Freq)
print TransferMatrix(ParList, Freq)
print TransferMatrix(ParList, Freq)[0,0]
print np.array([[1./TransferMatrix(ParList, Freq)[0, 0]], [0.] ])
ee = np.dot(TransferMatrix_After(ParList, Freq, 1.5e-3), np.array([[1./TransferMatrix(ParList, Freq)[0, 0]], [0.] ]) )
print ee[0] + ee[1]
'''


#plt.plot(DistList/1e-3, EFieldList.real)
#plt.plot(DistList/1e-3, EFieldList.imag)
#plt.plot(DistList/1e-3, np.sqrt(EFieldList.real**2 + EFieldList.imag**2), label = 'E')
#plt.plot(DistList/1e-3, BFieldList.real)
#plt.plot(DistList/1e-3, BFieldList.imag)
#plt.plot(DistList/1e-3, np.sqrt(BFieldList.real**2 + BFieldList.imag**2), label = 'B')
#plt.plot(DistList/1e-3, DFieldList.real)
#plt.plot(DistList/1e-3, DFieldList.imag)
#plt.plot(DistList/1e-3, np.sqrt(DFieldList.real**2 + DFieldList.imag**2), label = 'D')
#plt.plot(DistList/1e-3, HFieldList.real)
#plt.plot(DistList/1e-3, HFieldList.imag)
plt.plot(DistList/1e-3, np.sqrt(HFieldList.real**2 + HFieldList.imag**2), label = 'H')
plt.legend()
plt.savefig("FieldDist.pdf")