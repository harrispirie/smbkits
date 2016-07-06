import pylab as plt
import matplotlib as mpl
import numpy as np


color = ['black', 'violet', 'blue','cyan', 'green', 'yellow','gold', 'orange','red']

def plot(spec, temperatures = ['5K','10K'], greensFunctions = True, fitting = True, transport = True, harry = False):

    if greensFunctions:
        ix = 0
        fig = plt.figure(figsize = [10,10])
        for key in temperatures:
            plt.subplot(3,3,ix+1)
            plt.title(str(int(spec[key].T))+'K')
            plt.pcolormesh(np.linspace(0,2,1e4), spec[key].en, spec[key]._Gff+spec[key]._Gcc, cmap = mpl.cm.gist_heat_r)
            plt.ylim(spec[key].en[0],spec[key].en[-1])
            if ix != 6:
                plt.tick_params(labelbottom = 'off', labelleft = 'off')
            else:
                plt.xlabel('Momentum')
                plt.ylabel('Energy (meV)')
            ix += 1

    if fitting:
        ix = 0
        fig = plt.figure(figsize = [12,8])
        for key in temperatures:
            plt.subplot(2,3,1)
            plt.title('Hybridization')
            plt.plot(spec[key].T, spec[key].X[3], 'ob' , ms = 8, mec = 'k')
            #plt.ylim(112,114)
            plt.subplot(2,3,2)
            plt.title('f-electron self energy')
            plt.plot(spec[key].T, spec[key].X[1], 'ob', ms = 8, mec = 'k')
            plt.subplot(2,3,3)
            plt.title('Conductivity')
            plt.plot(spec[key].T,spec[key].conductivity[0], 'ob',  ms = 8, mec = 'k')
            plt.plot(spec[key].T, spec[key].thermalBroadening[0], 'or', ms = 8, mec = 'k')
            plt.subplot(2,3,4)
            plt.title('Tunneling ratio')
            plt.plot(spec[key].T,spec[key].X[4], 'ob', ms = 8, mec = 'k')
            plt.subplot(2,3,5)
            plt.title('Upper band minimum')
            plt.plot(spec[key].T,spec[key].upperBand, 'ob', ms = 8, mec = 'k')
            plt.subplot(2,3,6)
            plt.title('f-electron energy')
            plt.plot(spec[key].T,spec[key].X[2], 'ob', ms = 8, mec = 'k')
            ix += 1

    if transport:
        ix = 0
        fig = plt.figure(figsize = [12,8])
        for key in temperatures:
            plt.subplot(2,3,1)
            plt.title('Experiment')
            plt.plot(spec[key].en_crop,spec[key].didv_fit ,color = color[ix%len(color)], lw = 1)
            plt.plot(spec[key].en,spec[key].didv/spec[key].didv[-1] ,'ok', ms = 1)
            plt.subplot(2,3,2)
            plt.title('d-band')
            plt.plot(spec[key].en_crop,spec[key].fit_C,color = color[ix%len(color)], lw = 1)
            plt.subplot(2,3,3)
            plt.title('f-band')
            plt.plot(spec[key].en_crop,spec[key].fit_F/(spec[key].X[4]**2),color = color[ix%len(color)], lw = 1)
            plt.subplot(2,3,4)
            plt.title('Combined DOS')
            plt.plot(spec[key].xVector,spec[key].fnCombinedDOS(spec[key].xVector),color = color[ix%len(color)], lw = 1)
            plt.plot(-spec[key].kbT,spec[key].fnCombinedDOS(-spec[key].kbT),'o',color = color[ix%len(color)],ms = 5)
            plt.plot(spec[key].kbT,spec[key].fnCombinedDOS(spec[key].kbT),'o',color = color[ix%len(color)],ms = 5)
            plt.subplot(2,3,5)
            plt.title('Interference')
            plt.plot(spec[key].en_crop, spec[key].fit_CF / (2*spec[key].X[4]), color = color[ix%len(color)], lw = 1)
            plt.subplot(2,3,6)
            plt.title('Difference')
            plt.plot(spec[key].en, spec[key].difference, color = color[ix%len(color)], lw = 1)
            ix += 1

    if harry:
        ix = 0
        fig = plt.figure(figsize = [12,8])
        for key in temperatures:
            plt.subplot(2,3,1)
            plt.title('Bare f-band')
            plt.pcolormesh(np.linspace(0,2,1e4), spec[temperatures[0]].en, spec[temperatures[0]].Ekf)
            plt.ylim(-100, 100)
            plt.subplot(2,3,2)
            plt.title('Bare c-band')
            plt.pcolormesh(np.linspace(0,2,1e4), spec[temperatures[0]].en, spec[temperatures[0]].Ekc)
            plt.ylim(-100,100)
            plt.subplot(2,3,3)
            plt.title('Flat DOS')
            plt.plot(spec[key].en, spec[key].flatDidv, '-', color = color[ix%len(color)], lw=1)
            plt.xlim(-100,100)
            plt.subplot(2,3,4)
            plt.title('Gcc0')
            ix += 1
            

    plt.show()


