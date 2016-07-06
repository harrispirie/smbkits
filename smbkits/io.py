from scipy.interpolate import interp1d
import cPickle as pickle


def save(spec, filename):
	_removeFnCombinedDOS(spec)
	with open(filename, 'w') as output:
		pickle.dump(spec,output,-1)

def load(filename):
	with open(filename, 'r') as input:
		spec = pickle.load(input)
	_addFnCombinedDOS(spec)
	return spec


#  Need to remove function instance from class before pickling.
def _addFnCombinedDOS(spec):
        for key in spec.keys():
                spec[key].fnCombinedDOS = interp1d(spec[key].en_crop, spec[key].fit_C + spec[key].fit_F / (spec[key].X[4]**2), kind = 'cubic')

def _removeFnCombinedDOS(spec):
        for key in spec.keys():
                del spec[key].fnCombinedDOS

