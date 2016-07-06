import numpy as np
import scipy.io as sio
import nanoio as nio
import os

class Spectra(object):
    def __init__(self,didv,en):
        self.didv = didv
        self.en = en
        
class DOSmap(object):
    def __init__(self,LIY,en,Z):
        didvAvg = [np.mean(LIY[ix]) for ix,enVal in enumerate(en)]
        self.didv = didvAvg/didvAvg[-1]
        self.LIY = LIY/didvAvg[-1]
        self.en = en
        self.Z  = Z

dataFolder = '/Users/Harry/Documents/University/Hoffman Lab/SmB6/Data/_Good Data/'
validMatFiles = []; valid3dsFiles = []
for fileID in os.listdir(dataFolder):
    if fileID.endswith('.mat'): validMatFiles.append(fileID[:-4])
    if fileID.endswith('.3ds'): valid3dsFiles.append(fileID[:-4])

def importSpectra(fileName):
    ''' Import '.mat' file type as spectra object.  
    Automatically adjusts for bias offset by finding the zero crossing of the current.
    Normalizes the spectra by the end point.
    
    Usage: spec8K = importSpectra('spec8K')
    '''
    try: matFile = sio.loadmat( dataFolder + fileName + '.mat')
    except:
        print 'ERR: File not found. Must be one of the following files:'
        for fileID in validMatFiles: print fileID
        return 0
    try:
        matEn = matFile['data'][:,0]*1000
        matCurr = matFile['data'][:,1]
        matDidv = matFile['data'][:,2]
        biasOffset = matEn[np.argmin(abs(matCurr))]
        correctEn = matEn - biasOffset
        return Spectra(matDidv/matDidv[-1], correctEn)
    except:
        print 'File not in standard format, attempting to read...'
        try:
            matEn = matFile['bias'][:,0]
            matDidv = matFile['data'][:,0]
            print 'WARNING: No current channel found. Bias offset not determined.'
            return Spectra(matDidv/matDidv[-1],matEn)
        except: print 'ERR: Not standard keys. Resort to manual definitions with keys:\n{:}'.format(matFile.keys())


 

def importDOSmap(fileName):
    ''' Import '.3ds' file type as DOSmap object.
    Automatically adjusts for bias offset by finding the zero crossing of the average current.
    Normalizes the DOS map by the end point in the average didv.
    
    Usage: QPI_Gd = importDOSmap('QPI_Gd_001')
    '''
    try: grid = nio.Grid(dataFolder + fileName + '.3ds')
    except:
        print 'ERR: File not found. Must be one of the following files:'
        for fileID in valid3dsFiles: print fileID
        return 0
    try:
        avgCurrent = [np.mean(grid.I[ix]) for ix,en in enumerate(grid.en)]
        biasOffset = grid.en[np.argmin(np.abs(avgCurrent))]
        return DOSmap(grid.LIY,grid.en - biasOffset,grid.Z)
    except:
        print 'ERR: File not in standard format for processing.  Returned Grid object for manual processing.'
        return grid



    
    
    
    
    
    