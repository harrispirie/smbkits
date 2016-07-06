import numpy as np
import MorrFit
from scipy.optimize import minimize

class TemperatureSpectra(MorrFit.Didv):
    def __init__(self, didv, en):
        ###   NOTE: IM NORMALIZING DIDV FROM THE GET-GO   ###
        self.didv = didv/didv[-1]
        self.en = en
        self.N = len(en)
        self.res = 1e4
        
    def fitPolyBackground(self,fitRange,deg = 1):
        self.en_bk, self.didv_bk = self._HPmegaCrop(didv=self.didv,fitRange=fitRange)
        polyCoeff = np.polyfit(self.en_bk,self.didv_bk,deg)
        self.polyBackgroundFunction = np.poly1d(polyCoeff)
        self.fit_bk = self.polyBackgroundFunction(self.en)
        self.didvCorrected = self.didv - self.fit_bk
   
    def fitMorrSpec(self, didv, fitRange, X0 = None, method = 'SLSQP'):
        self.en_crop, self.didv_crop = self._HPmegaCrop(didv=didv, fitRange=fitRange)
        bounds = [(0,None),(0,None),(None,None),(0,None),(None,None),(None,None)]
        if X0 is None:
            X0 = np.array([1,2.83,-8.38,105.3,-0.026,0])
        result = minimize(self._HPchi4,X0,method = method,bounds = bounds)
        self.X = result.x
        
    def fitMorrDispersive(self, didv, fitRange, X0 = None, method = 'SLSQP'):
        self.en_crop, self.didv_crop = self._HPmegaCrop(didv=didv, fitRange=fitRange)
        bounds = [(0,None),(0,None),(-10,-1.9),(0,None),(None,None),(0,5),(None,None)]
        if X0 is None:
            X0 = np.array([1.95,6.65,-3.15,29.40,-0.006,6.72,-4.051])
        result = minimize(self._HPchiDispersive,X0,method = method,bounds = bounds)
        self.X = result.x
        
    def fitMorrDispersiveSmBFix(self, didv, fitRange, X0 = None, method = 'SLSQP'):
        self.en_crop, self.didv_crop = self._HPmegaCrop(didv=didv, fitRange=fitRange)
        bounds = [(0,None),(0,None),(-10,-1.9),(0,None),(None,None),(0,5),(-1601,-1599),(None,None)]
        if X0 is None:
            X0 = np.array([1.95,6.65,-3.15,29.40,-0.006,6.72,-1600,-4])
        result = minimize(self._HPchiDispersiveSmBFix,X0,method = method,bounds = bounds)
        self.X = result.x
    
    def fitMorrNorm(self, didv, fitRange, X0 = None, method = 'SLSQP'):
        self.en_crop, self.didv_crop = self._HPmegaCrop(didv=didv, fitRange=fitRange)
        bounds = [(0,None),(0,None),(None,None),(0,None),(None,None),(None,None),(None,None)]
        if X0 is None:
            X0 = np.array([1,2.83,-8.38,105.3,-0.026,0,0])
        result = minimize(self._HPchi5,X0,method = method,bounds = bounds)
        self.X = result.x

    def _HPmegaCrop(self, didv, fitRange):
        'fitRange is a list [start,stop,start,stop,...] which can delete sections of the data and the energy.'
        fitIndex = []
        enRange = []
        didvRange = []
        for energy in fitRange:
            a = [ix for ix in range(len(self.en)) if self.en[ix] >= energy]
            if a == []:
                if energy <= self.en[0]:
                    a = [0]
                else:
                    a = [len(self.en)-1]
            fitIndex.append(a[0])
        startIndex = [fitIndex[ix] for ix in range(len(fitIndex)) if ix % 2 is 0]
        stopIndex  = [fitIndex[ix] for ix in range(len(fitIndex)) if ix % 2 is 1]
        for start in startIndex:
            stop = stopIndex[startIndex.index(start)]
            enRange += self.en.tolist()[start:stop+1]
            didvRange += didv.tolist()[start:stop+1]
        return np.array(enRange), np.array(didvRange)
    
    def bandstructure(self, minEn, maxEn, N=1e3):
        en = np.linspace(minEn, maxEn, N)
        self._HPmodel(self.X[0],self.X[1],self.X[2],self.X[3],self.X[4],en)
        self.bands = np.imag(self._Gff) + np.imag(self._Gcc)
        self.combinedDOS = np.sum(self.Nc.T, axis=0) + np.sum(self.Nf.T, axis=0)
        
    
    def bandstructureDispersive(self, minEn, maxEn, N=1e3):
        en = np.linspace(minEn, maxEn, N)
        self.k = np.linspace(-2,2,self.res)
        y,c,f,cf = self._MorrDispersive(self.X[0],self.X[1],self.X[2],self.X[3],self.X[4],self.X[5],en)
        self.bands = np.imag(self.Gff) + np.imag(self.Gcc)
        self.enHR = en
        self.didvHR = y + self.X[6]
        self.combinedDOS = np.sum(self.Nc.T, axis=0) + np.sum(self.Nf.T, axis=0)
        self.Ef = self._Ef(self.k,self.X[5],self.X[2])
        self.Ec = self._Ec(self.k)
        self.Eup = self._Eup(self.k,self.X[5],self.X[2],self.X[3])
        self.Edw = self._Edw(self.k,self.X[5],self.X[2],self.X[3])

        
    def _MorrDispersive(self,rc,rf,ef,v,ratio,amplitude,en):
        N = self.res;
        r = np.linspace(0,2,N)
        r2 = (r**2.-1.)*1600.
        r3 = amplitude*np.cos(r*np.pi/2.0) + ef
        self.Ekc = np.tile (r2,[len(en),1])
        self.Ekf = np.tile (r3,[len(en),1])
        omega = np.tile(en,[N,1]).T
        Gcc0 = np.zeros([len(en),N])
        Gff0 = Gcc0; Gcc  = Gcc0; Gff  = Gcc0; Gcf  = Gcc0

        Gcc0 = (omega + rc*1j - self.Ekc)**(-1)
        Gff0 = (omega + rf*1j - self.Ekf)**(-1)
        self.Gcc  = (Gcc0**(-1) - v**2*Gff0)**(-1)
        self.Gff  = (Gff0**(-1) - v**2*Gcc0)**(-1)
        self.Gcf  = Gcc0*v*self.Gff;

        self.Nc = -np.imag(self.Gcc)*np.tile((r**2),[len(en),1])
        self.Nf = -np.imag(self.Gff)*np.tile((r**2),[len(en),1])
        self.Ncf= -np.imag(self.Gcf)*np.tile((r**2),[len(en),1])
        
        C = np.sum(self.Nc.T,axis = 0)
        F = ratio*ratio*np.sum(self.Nf.T,axis = 0)
        CF = 2*ratio*np.sum(self.Ncf.T,axis = 0)
        y = C + F + CF
        return y, C, F, CF
    
    def _HPSuperExtendDispersive(self, minEn, maxEn, N = 1e2):
        en = np.linspace(minEn, maxEn, N)
        f, self.superExtend_C, self.superExtend_F, self.superExtend_CF = self._MorrDispersive(self.X[0],self.X[1],self.X[2],
                                                                                       self.X[3],self.X[4],self.X[5],en)
        ft = f + self.X[6]
        self.didv_superExtend = ft
        self.en_superExtend = en
                
    def _HPchiDispersive(self,X):
        didv_fit, self.fit_C, self.fit_F, self.fit_CF = self._MorrDispersive(X[0],X[1],X[2],X[3],X[4],X[5],self.en_crop) 
        self.didv_fit = didv_fit + X[6]
        err = abs(self.didv_fit - self.didv_crop)
        self.chi2 = np.log(sum(err**2))
        return self.chi2
    
    def _HPchiDispersiveSmBFix(self,X):
        didv_fit, self.fit_C, self.fit_F, self.fit_CF = self._MorrDispersiveSmBFix(X[0],X[1],X[2],                                                                           X[3],X[4],X[5],X[6],self.en_crop) 
        self.didv_fit = didv_fit + X[7]
        err = abs(self.didv_fit - self.didv_crop)
        self.chi2 = np.log(sum(err**2))
        return self.chi2
    
    def _MorrDispersiveSmBFix(self,rc,rf,ef,v,ratio,amplitude,bottom,en):
        N = self.res;
        r = np.linspace(0,2,N)
        r2 = (r**2.-1.)*bottom
        r3 = amplitude*np.cos(r*np.pi/2.0) + ef
        self.Ekc = np.tile (r2,[len(en),1])
        self.Ekf = np.tile (r3,[len(en),1])
        omega = np.tile(en,[N,1]).T
        Gcc0 = np.zeros([len(en),N])
        Gff0 = Gcc0; Gcc  = Gcc0; Gff  = Gcc0; Gcf  = Gcc0

        Gcc0 = (omega + rc*1j - self.Ekc)**(-1)
        Gff0 = (omega + rf*1j - self.Ekf)**(-1)
        self.Gcc  = (Gcc0**(-1) - v**2*Gff0)**(-1)
        self.Gff  = (Gff0**(-1) - v**2*Gcc0)**(-1)
        self.Gcf  = Gcc0*v*self.Gff;

        self.Nc = -np.imag(self.Gcc)*np.tile((r**2),[len(en),1])
        self.Nf = -np.imag(self.Gff)*np.tile((r**2),[len(en),1])
        self.Ncf= -np.imag(self.Gcf)*np.tile((r**2),[len(en),1])
        
        C = np.sum(self.Nc.T,axis = 0)
        F = ratio*ratio*np.sum(self.Nf.T,axis = 0)
        CF = 2*ratio*np.sum(self.Ncf.T,axis = 0)
        y = C + F + CF
        return y, C, F, CF
    
    def _HPSuperExtend(self, minEn, maxEn, N = 1e2):
        en = np.linspace(minEn, maxEn, N)
        f, self.superExtend_C, self.superExtend_F, self.superExtend_CF = self._HPmodel(self.X[0],self.X[1],self.X[2],
                                                                                       self.X[3],self.X[4],en)
        ft = f + self.X[5]*np.linspace(self._HPBackgroundRange(minEn),self._HPBackgroundRange(maxEn),N) + self.X[6]
        #normIndex = [ix for ix in range(len(en)) if en[ix] > self.en[-1]]
        self.didv_superExtend = ft#/ft[normIndex[0]]
        self.en_superExtend = en
    
    def _HPSuperExtend2(self, minEn, maxEn, N = 1e2):
        en = np.linspace(minEn, maxEn, N)
        f, self.superExtend_C, self.superExtend_F, self.superExtend_CF = self._HPmodel(self.X[0],self.X[1],self.X[2],
                                                                                       self.X[3],self.X[4],en)
        ft = f + self.X[5]
        self.didv_superExtend = ft
        self.en_superExtend = en
    
    def _HPExtend(self):
        '''Extends a fit to span the energy range of the data set.'''
        f, self.superExtend_C, self.superExtend_F, self.superExtend_CF = self._HPmodel(self.X[0],self.X[1],self.X[2],
                                                                                       self.X[3],self.X[4],self.en)
        self.didv_extend = f + self.X[5]
      
    def _HPchi4(self,X):
        didv_fit, self.fit_C, self.fit_F, self.fit_CF = self._HPmodel(X[0],X[1],X[2],X[3],X[4],self.en_crop) 
        self.didv_fit = didv_fit + X[5]
        err = abs(self.didv_fit - self.didv_crop)
        self.chi4 = np.log(sum(err**2))
        return self.chi4
    
    def _HPchi5(self,X):
        didv_fit, self.fit_C, self.fit_F, self.fit_CF = self._HPmodel(X[0],X[1],X[2],X[3],X[4],self.en_crop) 
        self.didv_fit = didv_fit + X[5]*np.linspace(0,1,len(self.en_crop)) + X[6]
        err = abs(self.didv_fit - self.didv_crop)
        self.chi4 = np.log(sum(err**2))
        return self.chi4
    
    def _Ec(self,k):
        return (k**2.-1.)*1600.
    def _Ef(self,k,a,ef):
        return a*np.cos(k*np.pi/2.)+ef
    def _Eup(self,k,a,ef,v):
        return 0.5*(self._Ec(k) + self._Ef(k, a, ef)) + np.sqrt(0.25*(self._Ec(k) - self._Ef(k, a, ef))**2. + v**2.)
    def _Edw(self,k,a,ef,v):
        return 0.5*(self._Ec(k) + self._Ef(k, a, ef)) - np.sqrt(0.25*(self._Ec(k) - self._Ef(k, a, ef))**2. + v**2.)



