'''
The original VLE_model_fit has been modified to work specifically for Helium data
'''

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize  

font = {'size' : '24'}
plt.rc('font',**font)

Rg = 8.3144598e-5 #[m3 bar / K / mol]

class VLE(object):
    
    def __init__(self, Tsat, rhol, rhov, Psat, Mw=None, beta = 0.355):
        
        # Ensures that the lowest density IC corresponds to maximum temperature 
        rhol = rhol[np.argmax(Tsat):]
        rhov = rhov[np.argmax(Tsat):]
        Psat = Psat[np.argmax(Tsat):]
        Tsat = Tsat[np.argmax(Tsat):]  
    
        # Ensures that no Psats of 0
        Tsat = Tsat[Psat>0]
        rhol = rhol[Psat>0]
        rhov = rhov[Psat>0]
        Psat = Psat[Psat>0]
        
        self.Tsat = Tsat
        self.rhol = rhol
        self.rhov = rhov
        self.Psat = Psat
        self.logPsat = np.log10(Psat)
        self.invTsat = 1000./Tsat
        self.beta = beta
        self.Mw = Mw
        self.guessTc()
        self.fitRectScale()
        self.fitrhol()
        self.fitrhov()
        self.fitAntoine()
        self.Pc = self.PsatHat(self.Tc)
        self.Zc = self.computeZc(self.rhoc,self.Pc,self.Tc)
                
    def SSE(self,data,model):
        SE = (data - model)**2
        SSE = np.sum(SE)
        return SSE 

    def WSSE(self,data,model,weights):
        WSE = weights * (data - model)**2
        WSSE = np.sum(WSE)
        return WSSE

    def computew8reps(self,TsatAll,dataAll):
        w8All = np.zeros(len(TsatAll))
        TsatUnique = np.unique(TsatAll)
        for iTsat, Tsat in enumerate(TsatUnique):
            data = dataAll[TsatAll==Tsat]
            if len(data) > 1:
                w8 = (np.std(data))**2.
                w8All[TsatAll == Tsat] = w8      
        return w8All
    
    def computew8model(self,TsatAll,dataType,Tc,data=None):
        b012Dic = {'rhor':np.array([5.00e-4, 225e-14, 22.48]),'rhos':np.array([3.99e-4,2.24e-14,26.55]),'logPsat':np.array([3.51e-4,32.3e-14,24.86])}
        b012 = b012Dic[dataType]
        TrAll = TsatAll / Tc
        sigAll = b012[0] + b012[1]*np.exp(b012[2]*TrAll)
        if dataType == 'logPsat': 
            '''
            siglogPsat = 0.434 sigPsat / Psat
            sigPsat = RTsat sigrhov
            siglogPsat = 0.434 RTsat sigrhov / Psat
            siglogPsat = 0.434 sigrhov * rhov
            '''
            sigAll = 0.434 * sigAll * data #How to compute sig for f = log10(x)
        w8All = sigAll**2.
        return w8All
    
    def logPAntoine(self,b,T):
        logP = b[0] - b[1]/(b[2] + T)
        return logP

    def guessAntoine(self):
        slope, intercept, r_value, p_value, std_err = stats.linregress(self.invTsat,self.logPsat)
        guess = np.array([intercept,-slope*1000.,0])
        return guess
    
    def fitAntoine(self,w8Type='model'):
        Tfit = self.Tsat
        logPfit = self.logPsat
        
        if w8Type == 'uniform':
            w8logP = np.ones(len(Tfit))
        elif w8Type == 'replicates':
            w8logP = self.computew8reps(Tfit,logPfit)
        elif w8Type == 'model':
            w8logP = self.computew8model(Tfit,'logPsat',self.Tcguess,self.rhov)

        WSSElog = lambda b: self.WSSE(logPfit,self.logPAntoine(b,Tfit),w8logP)
        guess = self.guessAntoine()
        if len(Tfit) >= 3:
            bopt = minimize(WSSElog,guess).x
        else:
            bopt = guess
        self.boptPsat = bopt
        return bopt
    
    def logPsatHat(self,T):
        bopt = self.boptPsat
        logPsatHat = self.logPAntoine(bopt,T)
        return logPsatHat

    def PsatHat(self,T):
        PsatHat = 10.**(self.logPsatHat(T))
        return PsatHat
    
    def rholRectScale(self,b,T):
        beta = self.beta
        rhol = b[0] + b[1]*(b[2] - T) + b[3]*(b[2] - T)**beta
        return rhol
    
    def rhovRectScale(self,b,T):
        beta = self.beta
        rhov = b[0] + b[1]*(b[2] - T) - b[3]*(b[2] - T)**beta
        return rhov
    
    def rhorRect(self,b,T):
        rhor = b[0] + b[1]*(b[2]-T)
        return rhor
    
    def rhosScale(self,b,T,beta = 0):
        if beta == 0: beta = self.beta
        rhos = b[3]*(b[2]-T)**beta
        return rhos
    
    def guessRectScale(self):
        self.rhor = (self.rhol + self.rhov)/2.
        rhocGuess = np.mean(self.rhor)
        TcGuess = np.max(self.Tsat)/0.85
        guess = np.array([rhocGuess,2,TcGuess,50])
        ylin = self.rhol - rhocGuess #Modify to get decent guesses
        xlin = lambda b: self.rholRectScale(np.array([rhocGuess,b[0],TcGuess,b[1]]),self.Tsat) - rhocGuess
        SSEguess = lambda b: self.SSE(ylin,xlin(b))
        guess = np.array([2,50])
        bnd = ((0,None),(0,None))
        bopt = minimize(SSEguess,guess,bounds=bnd).x
        guess = np.array([rhocGuess,bopt[0],TcGuess,bopt[1]])
        return guess
    
    def guessTc(self):
        Tfit = self.Tsat
        rhorfit = (self.rhol + self.rhov)/2.
        rhosfit = (self.rhol - self.rhov)/2.
        w8rhor = self.computew8reps(Tfit,rhorfit)
        w8rhos = self.computew8reps(Tfit,rhosfit)
        WSSErhor = lambda b: self.WSSE(rhorfit,self.rhorRect(b,Tfit),w8rhor)
        WSSErhos = lambda b: self.WSSE(rhosfit,self.rhosScale(b,Tfit),w8rhos)
        WSSERectScale = lambda b: WSSErhor(b) + WSSErhos(b)
        guess = self.guessRectScale()
        if len(Tfit) >= 2:
            bnd = ((0,np.min(self.rhol)),(0,None),(np.max(Tfit),None),(0,None))
            bopt = minimize(WSSERectScale,guess,bounds=bnd).x
        else:
            bopt = guess
        self.Tcguess = bopt[2]
    
    def fitRectScale(self,w8Type='model'):
        Tfit = self.Tsat
        rhorfit = (self.rhol + self.rhov)/2.
        rhosfit = (self.rhol - self.rhov)/2.
    
        if w8Type == 'uniform':
            w8rhor = np.ones(len(Tfit))
            w8rhos = np.ones(len(Tfit))
        elif w8Type == 'replicates':
            w8rhor = self.computew8reps(Tfit,rhorfit)
            w8rhos = self.computew8reps(Tfit,rhosfit)
        elif w8Type == 'model':
            w8rhor = self.computew8model(Tfit,'rhor',self.Tcguess)
            w8rhos = self.computew8model(Tfit,'rhos',self.Tcguess)

        WSSErhor = lambda b: self.WSSE(rhorfit,self.rhorRect(b,Tfit),w8rhor)
        WSSErhos = lambda b: self.WSSE(rhosfit,self.rhosScale(b,Tfit),w8rhos)
        WSSERectScale = lambda b: WSSErhor(b) + WSSErhos(b)
        guess = self.guessRectScale()
        if len(Tfit) >= 2:
            bnd = ((0,np.min(self.rhol)),(0,None),(np.max(Tfit),None),(0,None))
            bopt = minimize(WSSERectScale,guess,bounds=bnd).x
        else:
            bopt = guess
        self.rhoc = bopt[0]
        self.Tc = bopt[2]
        self.boptRectScale = bopt
        return bopt
    
    def fitrhol(self):
        Tfit, rholfit = self.Tsat, self.rhol
        SSErhol = lambda b: self.SSE(rholfit,self.rholRectScale(b,Tfit))
        guess = self.boptRectScale
        #print(SSErhol(guess))
        if len(Tfit) >= 4: #rhol can get a better fit, although it extrapolates poorly
#            bnd = ((0,np.min(rholfit)),(0,None),(np.max(Tfit),None),(0,None))
            bnd = ((0,np.min(rholfit)),(0,None),(self.Tc,self.Tc),(0,None)) #Sets Tc to what RectScale got
            bopt = minimize(SSErhol,guess,bounds=bnd).x
        else:
            bopt = guess
        self.boptrhol = bopt
        #bopt = guess #If the optimization is not converging, this is a better option
        return bopt
    
    def rholHat(self,T):
        bopt = self.boptrhol
        rholHat = self.rholRectScale(bopt,T)
        return rholHat
    
    def fitrhov(self):
        Tfit, rhovfit = self.Tsat, self.rhov
        SSErhov = lambda b: self.SSE(rhovfit,self.rhovRectScale(b,Tfit))
        guess = self.boptRectScale
        #print(guess)
        #print(SSErhol(guess))
        if len(Tfit) >= 6: #rhov has problems fitting the data, better to just use from Rect Scale
            bnd = ((0,np.min(rhovfit)),(0,None),(np.max(Tfit),None),(0,None))
            bopt = minimize(SSErhov,guess,bounds=bnd).x
        else:
            bopt = guess
        self.boptrhov = bopt
        #bopt = guess #If the optimization is not converging, this is a better option
        return bopt
    
    def rhovHat(self,T):
        bopt = self.boptrhov
        #print(bopt)
        rhovHat = self.rhovRectScale(bopt,T)
        return rhovHat
    
    def rhorHat(self,T):
        bopt = self.boptRectScale
        rhorHat = self.rhorRect(bopt,T)
        return rhorHat    
    
    def rejectOutliers(self,data, m=2):
        return data[abs(data - np.mean(data)) < m * np.std(data)]
        
    def bootstrapCriticals(self,plothist=False,w8Type='model'):
        
        nBoots = 1000
        
        rhocBoots = np.zeros(nBoots)
        TcBoots = np.zeros(nBoots)
        PcBoots = np.zeros(nBoots)
        ZcBoots = np.zeros(nBoots)
        
        for iBoot in range(nBoots):
        
            ### First fit Tc and rhoc
            
            randbeta = np.random.uniform(0.325,0.385)
            
            randint = np.random.randint(0, len(self.Tsat),len(self.Tsat))
            Tfit = self.Tsat[randint]
            rholfit = self.rhol[randint]
            rhovfit = self.rhov[randint]
            rhorfit = (rholfit + rhovfit)/2.
            rhosfit = (rholfit - rhovfit)/2.
            logPfit = self.logPsat[randint]
                      
            if w8Type == 'uniform':
                w8rhor = np.ones(len(Tfit))
                w8rhos = np.ones(len(Tfit))
                w8logP = np.ones(len(Tfit))
            elif w8Type == 'replicates':
                w8rhor = self.computew8reps(Tfit,rhorfit)
                w8rhos = self.computew8reps(Tfit,rhosfit)
                w8logP = self.computew8reps(Tfit,logPfit)
            elif w8Type == 'model':
                w8rhor = self.computew8model(Tfit,'rhor',self.Tcguess)
                w8rhos = self.computew8model(Tfit,'rhos',self.Tcguess)
                w8logP = self.computew8model(Tfit,'logPsat',self.Tcguess,rhovfit)
                      
            WSSErhor = lambda b: self.WSSE(rhorfit,self.rhorRect(b,Tfit),w8rhor)
            WSSErhos = lambda b: self.WSSE(rhosfit,self.rhosScale(b,Tfit,beta=randbeta),w8rhos)
            WSSERectScale = lambda b: WSSErhor(b) + WSSErhos(b)
            guess = self.guessRectScale()
            if len(Tfit) >= 2:
                bnd = ((0,np.min(self.rhol)),(0,None),(np.max(Tfit),None),(0,None))
                bopt = minimize(WSSERectScale,guess,bounds=bnd).x
            else:
                bopt = guess
            rhocBoots[iBoot] = bopt[0]
            TcBoots[iBoot] = bopt[2]
            
            ### Then fit Pc
            
            WSSElog = lambda b: self.WSSE(logPfit,self.logPAntoine(b,Tfit),w8logP)
            guess = self.guessAntoine()
            if len(Tfit) >= 3:
                bopt = minimize(WSSElog,guess).x
            else:
                bopt = guess
                
            logPcBoot = self.logPAntoine(bopt,TcBoots[iBoot])
            PcBoot = 10.**(logPcBoot)
            
            PcBoots[iBoot] = PcBoot
                   
            ZcBoot = self.computeZc(rhocBoots[iBoot],PcBoots[iBoot],TcBoots[iBoot])
                   
            ZcBoots[iBoot] = ZcBoot
                   
        ### Removes outliers that might arise for poor bootstrapping
                   
        rhocBoots = self.rejectOutliers(rhocBoots)
        TcBoots = self.rejectOutliers(TcBoots)
        PcBoots = self.rejectOutliers(PcBoots)
        ZcBoots = self.rejectOutliers(ZcBoots)
        
        if plothist:
            plt.hist(rhocBoots,bins=50,color='k')
            plt.xlabel(r'$\rho_{\rm c}$ (kg/m$^3$)')
            plt.show()
    
            plt.hist(TcBoots,bins=50,color='k')
            plt.xlabel(r'$T_{\rm c}$ (K)')
            plt.show()
            
            plt.hist(PcBoots,bins=50,color='k')
            plt.xlabel(r'$P_{\rm c}$ (bar)')
            plt.show()
            
            plt.hist(ZcBoots,bins=50,color='k')
            plt.xlabel(r'$Z_{\rm c}$')
            plt.show()
                    
        urhoc = self.compute95CI(rhocBoots)
        uTc = self.compute95CI(TcBoots)
        uPc = self.compute95CI(PcBoots)
        uZc = self.compute95CI(ZcBoots)
                
        return urhoc, uTc, uPc, uZc
    
    def compute95CI(self,data):
        dataSorted = np.sort(data)
        i95 = int(0.025*len(dataSorted))
        dataTrimmed = dataSorted[i95:-i95] #Helps to avoid outliers also
        ### Assumes normal distribution of errors
        u95CI = 1.96 * np.std(dataTrimmed)
#        dataAvg = np.mean(data95)
#        dataError = (data95[-1] - data95[0])/(2*dataAvg)
        return u95CI
    
    def computeZc(self,rhoc,Pc,Tc):
        Vc = self.Mw / rhoc / 1000. #[m3/mol]
        Zc = Pc * Vc / Rg / Tc
        return Zc

def main():
    
    #### This is deprecated now that VLE requires passing the molecular weight to compute Zc
    
    # My Helium results:
    # Towhee
#    Tsat = np.array([6,7,8,9,10,11]) #[K]
#    rhol = np.array([0.31871202,0.30470999,0.289710787,0.273406827,0.254900377,0.231594877])*1000. #[kg/m3]
#    rhov = np.array([0.000179161,0.000713421,0.002167491,0.004683633,0.010768837,0.019240848])*1000. #[kg/m3]
#    Psat = np.array([2.221184667,10.17843633,34.23847833,80.274382,188.3929867,337.2498733])/100. #[bar]

    #Cassandra
    path = 'H:/Helium_ab_initio/Cassandra/Results/1400_10_md1/'
    VLCC = np.loadtxt(path+'VLCC.txt',skiprows=1)
    Psat = np.loadtxt(path+'Psat.txt',skiprows=1)[:,1]
    Tsat = VLCC[:,0]
    rhol = VLCC[:,1]
    rhov = VLCC[:,2]
    
##    #### Limit the range of data included in the fit
    TsatLow = 9
    TsatHigh = 15
    
    Psat = Psat[Tsat>=TsatLow]
    rhol = rhol[Tsat>=TsatLow]
    rhov = rhov[Tsat>=TsatLow]
    Tsat = Tsat[Tsat>=TsatLow]

    Psat = Psat[Tsat<=TsatHigh]
    rhol = rhol[Tsat<=TsatHigh]
    rhov = rhov[Tsat<=TsatHigh]
    Tsat = Tsat[Tsat<=TsatHigh]

    if Tsat[0] < Tsat[-1]:
        ###Have to make sure that Tsat[0] is the highest value since this code was written for ITIC
        Tsat = Tsat[::-1]
        rhol = rhol[::-1]
        rhov = rhov[::-1]
        Psat = Psat[::-1]
        
    #######

    He_fit = VLE(Tsat,rhol,rhov,Psat)
    
    invTsat = He_fit.invTsat/100 #For Helium it makes more sense to have 10/T
    logPsat = He_fit.logPsat
    Tsat = He_fit.Tsat
    rhol = He_fit.rhol
    rhov = He_fit.rhov
    rhor = He_fit.rhor
    Psat = He_fit.Psat
#    He_fit.fitRectScale()
#    #Tc = He_fit.boptRectScale[2] #Very poor Tc since only using rhol
    Tc = He_fit.Tc
    Pc = He_fit.Pc
    rhoc = He_fit.rhoc
    
#    
##    invTsat = 1000./Tsat
##    logPsat = np.log10(Psat)
#
    Tplot = np.linspace(min(Tsat),Tc,1000)
    invTplot = 10./Tplot
    logPsatplot = He_fit.logPsatHat(Tplot)
    rholplot = He_fit.rholHat(Tplot)
    rholRSplot = He_fit.rholRectScale(He_fit.boptRectScale,Tplot)
    Psatplot = He_fit.PsatHat(Tplot)
    rhovplot = He_fit.rhovHat(Tplot)
    rhovRSplot = He_fit.rhovRectScale(He_fit.boptRectScale,Tplot)
    rhorplot = He_fit.rhorHat(Tplot)
    Psatsmoothed = He_fit.PsatHat(Tsat)
    rholsmoothed = He_fit.rholHat(Tsat)
#    print(Psatplot)
#    print(Psatsmoothed)

    urhoc, uTc, uPc = He_fit.bootstrapCriticals(plothist=True)
    ulogPc = (np.log10(Pc+uPc) - np.log10(Pc-uPc))/2.
    uinvTc = (10./(Tc-uTc)-10./(Tc+uTc))/2.
    
    Mw_He = 4.0026 #[gm/mol]
             
    Pc_Kofke = 0.95*10. #[bar]
    rhoc_Kofke = 27.5 * Mw_He #[kg/m3]
    Tc_Kofke = 13.05 #[K]
    
    uPc_Kofke = 0.2*10. #[bar]
    urhoc_Kofke = 2.5 * Mw_He #[kg/m3]
    uTc_Kofke = 0.05 #[K]    

    ulogPc_Kofke = (np.log10(Pc_Kofke+uPc_Kofke) - np.log10(Pc_Kofke-uPc_Kofke))/2.
    uinvTc_Kofke = (10./(Tc_Kofke-uTc_Kofke)-10./(Tc_Kofke+uTc_Kofke))/2.     
    
    print('Critical temperature = '+str(np.round(Tc,3))+r'$ \pm$ '+str(np.round(uTc,3))+' K, Critical Pressure = '+str(np.round(Pc,3))+r'$ \pm$ '+str(np.round(uPc,3))+' bar, Critical Density = '+str(np.round(rhoc,1))+' $\pm$ '+str(np.round(urhoc,2))+' (kg/m$^3$).')
    
    plt.figure(figsize=[8,8])
    
    plt.plot(invTsat,logPsat-1,'ro',mfc='None',label='GEMC')
    plt.plot(invTplot,logPsatplot-1,'k',label='Fit, GEMC')
    plt.errorbar(10./Tc,np.log10(Pc)-1,xerr=uinvTc,yerr=ulogPc,fmt='b*',markersize=10,mfc='None',label='Critical, GEMC')
    plt.errorbar(10./Tc_Kofke,np.log10(Pc_Kofke)-1,xerr=uinvTc_Kofke,yerr=ulogPc_Kofke,mfc='None',fmt='gs',markersize=10,label='Critical, VEOS')
    plt.xlabel('10/T (K)')
    plt.ylabel(r'$\log_{10}$(P$_{\rm sat}$/MPa)')
    plt.legend()
    plt.show()
    
    plt.plot(Tsat,logPsat,'ro')
    plt.plot(Tplot,logPsatplot,'k')
    plt.errorbar(Tc,np.log10(Pc),xerr=uTc,yerr=ulogPc,fmt='b*')
    plt.xlabel('Temperature (K)')
    plt.ylabel('log(Psat/bar)')
    plt.show()
    
    plt.plot(Tsat,Psat,'ro')
#    plt.plot(Tsat,Psatsmoothed,'gx')
    plt.plot(Tplot,Psatplot,'k')
    plt.errorbar(Tc,Pc,xerr=uTc,yerr=uPc,fmt='b*')
    plt.errorbar(Tc_Kofke,Pc_Kofke,xerr=uTc_Kofke,yerr=uPc_Kofke,mfc='None',fmt='gs')
    plt.xlabel('Temperature (K)')
    plt.ylabel('Psat (bar)')
    plt.show()

    plt.figure(figsize=[8,8])

    plt.plot(rhol,Tsat,'ro',mfc='None',label='GEMC')
    plt.plot(rhov,Tsat,'ro',mfc='None')
    plt.plot(rhor,Tsat,'ro',mfc='None')
#    plt.plot(rholsmoothed,Tsat,'gx')
#    plt.plot(rholplot,Tplot,'g')
#    plt.plot(rhovplot,Tplot,'g')
    plt.plot(rhorplot,Tplot,'k',label='Fit, GEMC')
    plt.plot(rholRSplot,Tplot,'k')
    plt.plot(rhovRSplot,Tplot,'k')
    plt.errorbar(rhoc,Tc,xerr=urhoc,yerr=uTc,fmt='b*',markersize=10,mfc='None',label='Critical, GEMC')
    plt.errorbar(rhoc_Kofke,Tc_Kofke,xerr=urhoc_Kofke,yerr=uTc_Kofke,fmt='gs',markersize=10,mfc='None',label='Critical, VEOS')
    plt.xlabel(r'Density (kg/m$^3$)')
    plt.ylabel('Temperature (K)')
    plt.legend()
    plt.show()
    
    Rg = 8.3144598e-5 #[m3 bar / K / mol]
    Vv = Mw_He/rhov/1000. #[m3/mol]
    Zv = Psat*Vv/Rg/Tsat
    
    Zvplot = Psatplot * Mw_He / rhovRSplot / Rg / Tplot / 1000.
    
    Vc = Mw_He/rhoc/1000.
    Zc = Pc * Vc / Rg / Tc
    
    uZc = Zc * np.sqrt((uPc/Pc)**2. + (urhoc/rhoc)**2. + (uTc/Tc)**2.)
    
    Vc_Kofke = Mw_He/rhoc_Kofke/1000.
    Zc_Kofke = Pc_Kofke * Vc_Kofke / Rg / Tc_Kofke
    
    uZc_Kofke = Zc_Kofke * np.sqrt((uPc_Kofke/Pc_Kofke)**2. + (urhoc_Kofke/rhoc_Kofke)**2. + (uTc_Kofke/Tc_Kofke)**2.)
    
#    print(Zv,Zc,uZc)
    
    Tplot_Zv = Tplot[Zvplot>0]
    Zvplot = Zvplot[Zvplot>0]
    
    Tplot_Zv = Tplot_Zv[Zvplot<1]
    Zvplot = Zvplot[Zvplot<1]
    
    plt.figure(figsize=[8,8])

    plt.plot(Tsat,Zv,'ro',mfc='None',label='GEMC')
#    plt.plot(Tplot_Zv,Zvplot,'k',label='Fit, GEMC')
    plt.errorbar(Tc,Zc,xerr=uTc,yerr=uZc,fmt='b*',markersize=10,mfc='None',label='Critical, GEMC')
#    plt.errorbar(Tc_Kofke,Zc_Kofke,xerr=uTc_Kofke,yerr=uZc_Kofke,fmt='gs',mfc='None',label='Critical, Kofke')
    plt.ylabel(r'Vapor Compressibility Factor')
    plt.xlabel('Temperature (K)')
    plt.legend()
    plt.show()

if __name__ == '__main__':
    
    main()