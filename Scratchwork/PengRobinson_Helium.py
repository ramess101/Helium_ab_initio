# -*- coding: utf-8 -*-
"""
Created on Fri Nov 17 09:20:52 2017

@author: ram9
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import minimize, root

Rg = 8.3144598e-6 #[m3Mpa/K/mol]
Mw = 4.0026/1000 #[kg/mol]

Tc_PR = 11.73 #[K]
Pc_PR = 0.568 #[MPa]
omega_PR = 0.

class PengRob(object):
    
    def __init__(self, Tc, Pc, omega, Mw, Temp):
        self.Tc = Tc
        self.Pc = Pc
        self.omega = omega
        self.Mw = Mw
        self.Temp = Temp
        self.a = self.calc_a()
        self.b = self.calc_b()
        self.kappa = self.calc_kappa()
            
    def calc_a(self):
        """ Return the a constant """
        a = 0.45724*Rg**2*self.Tc**2/self.Pc
        return a

    def calc_b(self):
        """ Return the b constant """
        b = 0.07780 * Rg * self.Tc / self.Pc
        return b
    
    def calc_kappa(self):
        """ Return the kappa constant """
        kappa = 0.37464 + 1.54226 * self.omega - 0.26992 * self.omega**2
        return kappa
    
    def calc_alpha(self,Temp):
        """ Return the alpha temperature dependent constant """
        alpha = (1+self.kappa*(1-np.sqrt(Temp/self.Tc)))**2
        return alpha
    
    def calc_A(self,Temp,Press):
        """ Return the A temperature dependent constant """
        alpha = self.calc_alpha(Temp)
        A = self.a * alpha * Press / Rg**2 / Temp**2
        return A

    def calc_B(self,Temp,Press):
        """ Return the B temperature dependent constant """
        B = self.b * Press / Rg / Temp
        return B
    
    def calc_Z(self,Temp,Press,A,B):

        PR_Z = lambda Z: Z**3 - (1.-B)* Z**2. + (A-3.*B**2.-2.*B)*Z - (A*B - B**2. - B**3.)
        
        Zguess = 1
        opt = root(PR_Z,Zguess)
        Zv = opt.x[0]
        
        Zguess = 0.
        opt = root(PR_Z,Zguess)
        Zl = opt.x[0]
        
        return Zv, Zl
    
    def calc_f(self,Temp,Press):
        
        A = self.calc_A(Temp,Press)
        B = self.calc_B(Temp,Press)
        
        Zv, Zl = self.calc_Z(Temp,Press,A,B)
        
        fug = lambda Z: Press * np.exp((Z-1) - np.log(Z-B) - (A/(2*np.sqrt(2)*B))*np.log((Z+(1+np.sqrt(2))*B)/(Z+(1-np.sqrt(2))*B)))
               
        fv = fug(Zv)
        fl = fug(Zl)
        
        return fv, fl, Zv, Zl
    
    def calc_rho(self,Z,Temp,Press):

        rho = Press / Temp / Z / Rg * self.Mw #[kg/m3]
        
        return rho
    
    def calc_VLE(self,Temp_array):
        
        Psat = np.zeros(len(Temp_array))
        rhov = np.zeros(len(Temp_array))
        rhol = np.zeros(len(Temp_array))
        
        for iTemp, Temp in enumerate(Temp_array):
        
            Press = 0.1*self.Pc
            
            conv = False
            max_it = 100
            it = 0.
            
            while conv == False and it < max_it:
            
                fv,fl,Zv,Zl = self.calc_f(Temp,Press)
            
                if np.abs((fl/fv)-1.)< 1e-4:
                
                    conv = True
                    
                else:
                    
                    Press *= fl/fv
                    it += 1
                    
            Psat[iTemp] = Press
            rhov[iTemp] = self.calc_rho(Zv,Temp,Press)
            rhol[iTemp] = self.calc_rho(Zl,Temp,Press)
                
        return Psat, rhov, rhol
            
    
HePR = PengRob(Tc_PR,Pc_PR,omega_PR,Mw,np.array([Tc_PR,10]))

Tplot = np.linspace(6,0.99*Tc_PR,200)
Psatplot, rhovplot, rholplot = HePR.calc_VLE(Tplot)

plt.plot(rholplot,Tplot,'k')
plt.plot(rhovplot,Tplot,'k')
plt.show()

plt.plot(1000./Tplot,np.log10(Psatplot),'k')
plt.show()