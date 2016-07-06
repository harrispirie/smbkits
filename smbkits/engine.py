import numpy as np
import pylab as plt
import matplotlib as mpl
import scipy.io as sio
import scipy.constants as const
from scipy.interpolate import interp1d
from scipy.integrate import simps, quad
import smbfits
import MorrFit
import time

dir = '/Users/Harry/Documents/University/Hoffman Lab/SmB6/Data/Temperature Data/'

str0 = ['5K','10K']
str1 = ['5K','10K','15K','30K','50K']
str2 = ['2.2K', '4.4K', '8K', '13K', '20K', '25K', '30K1']
str3 = ['8Kg','15Kg','20Kg','25Kg','30Kg','35Kg','40Kg','45Kg','50Kg']

def processData(strix,start = -25, fitRange = [-25,100], method = 'SLSQP', run = True):
    ix = -1; t1 = time.clock()
    X0 = np.array([ 1.76120251e+00, 5.31298602e+00, -3.36660978e+00,
                 4.25126930e+01, -4.08585819e-03, -4.03481375e+00])
    X0 = np.array([ 1.94594665e+00, 6.65650865e+00, -6.00000000e+00, 
                    2.93997298e+01, -6.04674289e-03, 4, 0])
    spec = {}

    for key in strix:
        ix += 1
        matFile = sio.loadmat(dir + 'spec' + key + '.mat')
        matEn = matFile['data'][:,0]*1000
        matCurr = matFile['data'][:,1]
        biasOffset = matEn[np.argmin(abs(matCurr))]
        correctEn = matEn - biasOffset
        matDidv = matFile['data'][:,2]
        spec[key] = smbfits.TemperatureSpectra(matDidv, correctEn)
        if run:
            #spec[key].fitPolyBackground([-200,-50,50,200])
            # Change: fitMorrSpec --> fitMorrNorm, didv --> didvCorrected, _HPSuperExtend2 --> _HPSuperExtend
            #         and X0 --> append 0 (above) for background subtraction.
            #          fitMorrSpec --> fitMorrDispersive, _HPSuperExtend2 --> _HPSuperExtendDispersive, change X0
            
            spec[key].fitMorrDispersive(didv = spec[key].didv, fitRange = fitRange, method = method,X0 = None)
            spec[key]._HPSuperExtendDispersive(-100,100) 
            spec[key].upperBand = upperBandMinimum(spec[key].X[2],spec[key].X[3])
            spec[key].bandstructureDispersive(-100,100)
            X0 = spec[key].X

            #Restivity
            spec[key].T = float(key.split('K')[0])
            spec[key].kbT = spec[key].T * const.k / const.e * 1000 #in meV
            spec[key].fnCombinedDOS = interp1d(spec[key].en_crop, 
                                               spec[key].fit_C + spec[key].fit_F / (spec[key].X[4]**2), kind = 'cubic')
            spec[key].xVector = np.linspace(spec[key].en_crop[0],spec[key].en_crop[-1],500)
            spec[key].conductivity = quad(spec[key].fnCombinedDOS, -spec[key].kbT , spec[key].kbT)
        
    t2 = time.clock()
    minutes = int(np.floor((t2-t1)/60)); seconds = round(t2-t1-60*minutes,1)
    print 'Time -  ', minutes, ':', seconds
    return spec

def upperBandMinimum(ef,v):
    return -800+ef/2. + np.sqrt( (-800-ef/2.)**2 + v**2 )



