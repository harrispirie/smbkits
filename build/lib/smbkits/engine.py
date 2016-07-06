import numpy as np
import pylab as plt
import matplotlib as mpl
import scipy.io as sio
import scipy.constants as const
from scipy.interpolate import interp1d
from scipy.integrate import simps, quad
import fit
import time

dir = '/Users/Harry/Documents/University/Hoffman Lab/SmB6/Data/Temperature Data/'

str0 = ['5K','10K']
str1 = ['5K','10K','15K','30K','50K']	
str2 = ['2.2K', '4.4K', '8K', '13K', '20K', '25K', '30K1']
str3 = ['8Kg','15Kg','20Kg','25Kg','30Kg','35Kg','40Kg','45Kg','50Kg']

def processData(strix,start = -25, method = 'SLSQP'):
	ix = -1; t1 = time.clock()
	X0 = np.array([  1.00000000e+00,   2.83000000e+00,  -5.08424882e+00, 1.13346446e+02,  -2.10958055e-02,   3.68803982e-02, 1.39194387e+00])

	spec = {}

	for key in strix:
		ix += 1
		matFile = sio.loadmat(dir + 'spec' + key + '.mat')
		matEn = matFile['data'][:,0]*1000
		matCurr = matFile['data'][:,1]
		biasOffset = matEn[np.argmin(abs(matCurr))]
		correctEn = matEn - biasOffset
		matDidv = matFile['data'][:,2]
		spec[key] = fit.Didv(matDidv, correctEn)
		spec[key].fit_didv(start = start, X0 = X0, method = method)
		spec[key].adjustedDidv = spec[key].didv/spec[key].didv[-1] - spec[key].X[5]*np.linspace(0,spec[key].en[-1],len(spec[key].en)) - spec[key].X[6]
		spec[key]._HPExtend(spec[key].X,spec[key].en) 
		spec[key].difference = spec[key].didv_extend - spec[key].didv/spec[key].didv[-1] 
		spec[key].upperBand = upperBandMinimum(spec[key].X[3],spec[key].X[4])
#		X0 = spec[key].X

		# Restivity
		spec[key].T = float(key.split('K')[0])
		spec[key].kbT = spec[key].T * const.k / const.e * 1000 	#in meV
		spec[key].fnCombinedDOS = interp1d(spec[key].en_crop, spec[key].fit_C + spec[key].fit_F / (spec[key].X[4]**2), kind = 'cubic')
		spec[key].xVector = np.linspace(spec[key].en_crop[0],spec[key].en_crop[-1],500)
		spec[key].conductivity = quad(spec[key].fnCombinedDOS, -spec[key].kbT , spec[key].kbT)
		spec[key].thermalBroadening = quad(spec[strix[0]].fnCombinedDOS, -spec[key].kbT, spec[key].kbT)

		# Remove slope from things
		bkg = np.linspace(0,len(spec[key].en_crop),len(spec[key].en))*spec[key].X[5]+spec[key].X[6]
		spec[key].slope = bkg/bkg[-1]
		spec[key].flatDidv = spec[key].didv/spec[key].didv[-1] - spec[key].slope

		
	
	t2 = time.clock()
	minutes = int(np.floor((t2-t1)/60))
	seconds = round(t2-t1-60*minutes,1)
	print 'Time -  ', minutes, ':', seconds
	return spec	

def upperBandMinimum(v,ef):
	return -800+ef/2. + np.sqrt( (-800-ef/2.)**2 + v**2 )


