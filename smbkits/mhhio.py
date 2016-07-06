import scipy.io as sio

def load(filePath):
	fileObject = sio.loadmat(filePath)
	typeID = []
	dataObjects = {}
	for dataKey in fileObject.keys():
		if dataKey not in ['__header__','__version__','__globals__']:
			typeID.append(dataKey)
	for ID in typeID:
		dataObjects[ID] = pyData(fileObject,ID)
	return dataObjects

def savem(filePath,pyObject):
	pySave = {}
	for key in pyObject.keys():
		pyDict = vars(pyObject[key])
		pyDict['info'] = vars(pyObject[key].info)
		pySave[key] = pyDict
	sio.savemat(filePath, pySave)

class pyData(object):
	def __init__(self,fileObject,typeID):
		contents = fileObject[typeID]
		for name in contents.dtype.names:
			setattr(self, name, contents[name][0,0])
		self.info = pyInfo(self.info)

class pyInfo(object):
	def __init__(self,info) :
		for name in info.dtype.names:
			 setattr(self, name, info[name][0,0])

