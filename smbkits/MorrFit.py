import numpy as np
import time
from scipy.optimize import minimize


class Didv(object):
    def __init__(self,didv,en):
        self.didv = didv
        self.en = en
        self.N = len(en)

    def fit_didv(self, start = None, stop = None, method = 'SLSQP', X0 = None, maxiter = None, output = False):
        """ Crops the average didv to the range specified by start and stop and fits this to the Figgins
        Model by constrained minimisation.  Minimization methods are SLSQP (default), TNC or L-BFGS-B"""
        t1 = time.clock()
        if start is None:
            start = self.en[0]
        if stop is None:
            stop = self.en[-1]
        self._HPcrop(start = start, stop = stop)
        bounds = [(0,None),(0,None),(None,None),(0,None),(None,None),(None,None),(None,None)]
        if X0 is None:
            X0 = np.array([1,2.83,-8.38,105.3,-0.026,-0.01,10.08])
        self.X = X0
        if maxiter is not None:
            options = {'maxiter':maxiter}
        if maxiter is None:
            options = {}
        result = minimize(self._HPchi2,self.X,method = method,bounds = bounds, options = options)
        self.X = result.x
        c = self._HPchi2(result.x)
        t2 = time.clock()
        minutes = int(np.floor((t2-t1)/60))
        seconds = round(t2-t1-60*minutes,1)
        self.summary =  '''\nElasped time: {minutes:5.0f} minutes {seconds:2.0f} seconds \n\n\nCalculated Parameter Values:\n
 c-electron self energy:\t {rc:5.2f}\n f-electron self energy:\t {rf:5.2f} \n f-electron energy:\t    {ef:5.2f}
 Hybridization:\t\t     {v:5.2f}\n Hopping amplitude ratio:\t{t:5.2f} \n Background slope:\t     {a:5.2f}
 Background Offset:\t    {b:5.2f}\n Final chi value:\t{c:5.2f} \n {message:20s}\n
                        '''.format(minutes=minutes, seconds=seconds, rc=self.X[0], rf=self.X[1], ef=self.X[2], v=self.X[3], t=self.X[4], a=self.X[5], b=self.X[6], c=c, message=result.message)
        if output is True:
            print self.summary

    def _HPcrop(self, start = None, stop = None):
        """ Crops didv to the energy range indicated by start and stop.  
            Normalizes at the end point and creates new attributes self.en_crop and self.didv_crop."""
        self.en_crop = [self.en[ix] for ix in range(self.N) if self.en[ix]>=start and self.en[ix] <=stop]
        self.didv_crop = [self.didv[ix] / self.didv[-1] for ix in range(self.N) if self.en[ix]>=start and self.en[ix] <=stop]
       
    def _HPchi2(self,X):
        f, self.fit_C, self.fit_F, self.fit_CF = self._HPmodel(X[0],X[1],X[2],X[3],X[4],self.en_crop) 
        ft = f + X[5]*np.linspace(0,1,len(self.en_crop)) + X[6]
        self._HPNormConstant = ft[-1]
        self.didv_fit = ft / ft[-1]
        err = abs(self.didv_fit - self.didv_crop)
        self.chi2 = np.log(sum(err**2))
        return self.chi2

    def _HPmodel(self,rc,rf,ef,v,ratio,en):
        # Size of grid for calculation - resolution in k space.
        N = 1e4;
   
        # Calculate bandstructre from ARPES, assume spherical d band, non-dispersive f.
        # Take r to run from 0.0 to 2.0,  -- the matlab code runs from 2/N to 2.0
        r = np.linspace(0,2,N)
        r2 = (r**2-1)*1600
        self.Ekc = np.tile (r2,[len(en),1])
        self.Ekf = np.zeros([len(en),N])+ef
        omega = np.tile(en,[N,1]).T

        # Initialization
        Gcc0 = np.zeros([len(en),N])
        Gff0 = Gcc0
        Gcc  = Gcc0
        Gff  = Gcc0
        Gcf  = Gcc0

        # Calculation
        Gcc0 = (omega + rc*1j - self.Ekc)**(-1)
        Gff0 = (omega + rf*1j - self.Ekf)**(-1)
        Gcc  = (Gcc0**(-1) - v**2*Gff0)**(-1)
        Gff  = (Gff0**(-1) - v**2*Gcc0)**(-1)
        Gcf  = Gcc0*v*Gff;

        # Properties for debugging
        self._Gcc0  = Gcc0
        self._Gff0  = Gff0
        self._Gcc   = Gcc
        self._Gff   = Gff
        self._Gcf   = Gcf
        self._r     = r
        self._r2    = r2
        self._omega = omega
        self._en    = en
        self._N     = N

        self.Nc = -np.imag(Gcc)*np.tile((r**2),[len(en),1])
        self.Nf = -np.imag(Gff)*np.tile((r**2),[len(en),1])
        self.Ncf= -np.imag(Gcf)*np.tile((r**2),[len(en),1])
        
        C = np.sum(self.Nc.T,axis = 0)
        F = ratio*ratio*np.sum(self.Nf.T,axis = 0)
        CF = 2*ratio*np.sum(self.Ncf.T,axis = 0)
        y = C + F + CF
        return y, C, F, CF

    def _HPExtend(self, en):
        f, self.extend_C, self.extend_F, self.extend_CF = self._HPmodel(self.X[0],self.X[1],self.X[2],self.X[3],self.X[4],en)
        ft = f + self.X[5]*np.linspace(self._HPBackgroundRange(self.en[0]),1,len(en)) + self.X[6]
        self.didv_extend = ft / ft[-1]
        
    def _HPBackgroundRange(self,x):
        y = (x - self.en_crop[0]) / (self.en_crop[-1] - self.en_crop[0])
        return y
        
    def _HPSuperExtend(self, minEn, maxEn, N = 1e2):
        en = np.linspace(minEn, maxEn, N)
        f, self.superExtend_C, self.superExtend_F, self.superExtend_CF = self._HPmodel(self.X[0],self.X[1],self.X[2],self.X[3],self.X[4],en)
        ft = f + self.X[5]*np.linspace(self._HPBackgroundRange(minEn),self._HPBackgroundRange(maxEn),N) + self.X[6]
        normIndex = [ix for ix in range(len(en)) if en[ix] > self.en[-1]]
        self.didv_superExtend = ft/ft[normIndex[0]]
        self.en_superExtend = en
        
 
