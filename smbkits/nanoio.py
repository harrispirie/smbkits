from struct import unpack
import numpy as np

class Grid(object):
    def __init__(self,filePath):
        if self._load3ds(filePath):
            try: self.LIY = self.data['LIY 1 omega (A)']
            except (KeyError):
                print 'Does not have channel called LIY 1 omega (A).  Looking for average channel instead...'
                try: self.LIY = self.data['LIY 1 omega [AVG] (A)']; print 'Found it!'
                except (KeyError):
                    print 'ERR: Average channel not found, resort to manual definitions.  Found channels:\n {:}'.format(self.data.keys())
            try: self.I   = self.data['Current (A)']
            except (KeyError): self.I = self.data['Current [AVG] (A)']
            self.Z   = self.parameters['Scan:Z (m)']
        else: raise NameError('NoSuchFile')

    def _load3ds(self,filePath):
        try: fileObj = open(filePath,'rb')
        except: return 0
        self.header={}
        while True:
            line = fileObj.readline().strip()
            if line == ':HEADER_END:': break
            splitLine = line.split('=')
            self.header[splitLine[0]] = splitLine[1]

        self.info={	'params'	:	int(self.header['# Parameters (4 byte)']),
                    'paramName'	:	self.header['Fixed parameters'][1:-1].split(';') +
                                    self.header['Experiment parameters'][1:-1].split(';'),
                    'channels'	:	self.header['Channels'][1:-1].split(';'),
                    'points'	:	int(self.header['Points']),
                    'sizex'		:	int(self.header['Grid dim'][1:-1].split(' x ')[0]),
                    'sizey'		:	int(self.header['Grid dim'][1:-1].split(' x ')[0]),
                    'dataStart'	:	fileObj.tell()
                    }

        self.en = np.linspace(float(self.header['Bias Spectroscopy>Sweep Start (V)']),
                              float(self.header['Bias Spectroscopy>Sweep End (V)']),
                              int(self.header['Bias Spectroscopy>Num Pixel'])
                              ) * 1000
                             
        self.data = {}; self.parameters = {}
        for channel in self.info['channels']:
            self.data[channel] = np.zeros([self.info['points'],self.info['sizex'],self.info['sizey']])
        for name in self.info['paramName']:
            self.parameters[name] = np.zeros([self.info['sizex'],self.info['sizey']])

        for ix in range(self.info['sizex']):
            for iy in range(self.info['sizey']):
                for name in self.info['paramName']:
                    value = unpack('>f',fileObj.read(4))[0]
                    self.parameters[name][ix,iy] = value
                
                for channel in self.info['channels']:
                    for ie in range(self.info['points']):
                        value = unpack('>f',fileObj.read(4))[0]
                        self.data[channel][ie,ix,iy] =value

        dataRead = fileObj.tell()
        fileObj.read()
        allData = fileObj.tell()
        if dataRead == allData: print 'File import successful.'
        else: print 'ERR: Did not reach end of file.'
        fileObj.close()
        return 1








