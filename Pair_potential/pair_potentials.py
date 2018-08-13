# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:31:30 2017

@author: ram9
"""

import numpy as np
import matplotlib.pyplot as plt
import math
from scipy.misc import derivative
from scipy.optimize import minimize, basinhopping

hartree_to_K = 315775.13
bohr_to_A = 0.529177 

def TangToennies(n,x):
    sum_in = 0.
    for i in range(n+1):
        sum_in += x**i/math.factorial(i)
    fn = 1-np.exp(-x)*sum_in
    return fn

def long_range(R,Cn,N0,N1,zeta):
    sum_N0N1 = 0.
    n_range = np.linspace(N0,N1,num=N1-N0+1)
    for i,n in enumerate(n_range):
        fn = TangToennies(int(n),zeta*R)
        sum_N0N1 += fn*Cn[i]/R**n
    return sum_N0N1

def short_range(R,ak,I0,I1,Pik):
    i_range = np.linspace(I0,I1,num=I1-I0+1)
    sum_kM = 0.
    for k, ak_k in enumerate(ak): #Since k only serves as an index in equation, I resorted back to index 0 python syntax
        sum_I0I1 = 0.
        for ik, i in enumerate(i_range):
            sum_I0I1 += Pik[k,ik]*R**i
            #print(ik,k,Pik[k,ik],Pik[k,ik]*R**i)
        sum_kM += np.exp(-ak_k*R)*sum_I0I1
    return sum_kM

def V_general(R,Cn,N0,N1,zeta,ak,I0,I1,Pik):
    R_bohr = R/bohr_to_A #convert Angstrom to bohr
    V_h = short_range(R_bohr,ak,I0,I1,Pik) - long_range(R_bohr,Cn,N0,N1,zeta) # Hartree
    V_K = V_h * hartree_to_K
    return V_K

Cn_BO = np.array([1.460977838,0,14.11785737,0,183.691075,-76.72571,3372,-3808.3254,85340,-170700,2860000])
zeta_BO = 4.45565354
N0_BO = 6
N1_BO = 16

ak_BO = np.array([2.190769009,3.915583621,8.873143625])
Pik_BO = np.array([26.91509667,17.40159689,-3.265662082,0.146700434,-262.457651,457.3005084,666.6249475,-101.5411647,239.5425543,673.2233558,862.6518076,786.2126675]).reshape(3,4)
I0_BO = -1
I1_BO = 2

Cn_AD = np.array([0.0011445,0,0.006519,0,0.0668])
zeta_AD = 4.809717508
N0_AD = 6
N1_AD = 10

ak_AD = np.array([2.75067514,5.360852666,7.137312082])
Pik_AD = np.array([0.016003492,-0.005179152,0.003200214,-0.502322895,-0.05281783,0.29399979,-0.078022274,1.83457621,0.946016994]).reshape(3,3)
I0_AD = 0
I1_AD = 2

Cn_REL = np.array([-0.000035322,0,-0.0003434,0,-0.00355798])
zeta_REL = 5.002125476
N0_REL = 4
N1_REL = 8

ak_REL = np.array([2.554427875,4.268222066,6.318015453])
Pik_REL = np.array([0.011625994,-0.007008058,0.000783673,-0.573711011,0.397167158,-0.110854902,0.622980694,0.487089036,0.681223618,]).reshape(3,3)
I0_REL = 0
I1_REL = 2

Cn_QED = np.array([5.77235E-07,0,1.37784E-06,7.61188E-05])
zeta_QED = 5.437355143
N0_QED = 3
N1_QED = 6

ak_QED = np.array([2.713859341,4.85539105])
Pik_QED = np.array([0.001349522,-0.001023625,0.000356862,-0.000196866,-0.000961842,-0.001080147]).reshape(2,3)
I0_QED = 0
I1_QED = 2

assert np.abs(long_range(3.,Cn_BO,N0_BO,N1_BO,zeta_BO) - 0.024620752) < 1e-5, "Error in long_range calculator"
assert np.abs(short_range(3.,ak_BO,I0_BO,I1_BO,Pik_BO) - 0.036552452) < 1e-5, "Error in short_range calculator"
assert np.abs(V_general(3.*bohr_to_A,Cn_BO,N0_BO,N1_BO,zeta_BO,ak_BO,I0_BO,I1_BO,Pik_BO) - 3767.7341) < 1e-3, "Error in total calculator"
      
def V_BO_hat(R):
    V = V_general(R,Cn_BO,N0_BO,N1_BO,zeta_BO,ak_BO,I0_BO,I1_BO,Pik_BO)
    return V

def V_AD_hat(R):
    V = V_general(R,Cn_AD,N0_AD,N1_AD,zeta_AD,ak_AD,I0_AD,I1_AD,Pik_AD)
    return V

def V_REL_hat(R):
    V = V_general(R,Cn_REL,N0_REL,N1_REL,zeta_REL,ak_REL,I0_REL,I1_REL,Pik_REL)
    return V

def V_QED_hat(R):
    V = V_general(R,Cn_QED,N0_QED,N1_QED,zeta_QED,ak_QED,I0_QED,I1_QED,Pik_QED)
    return V  

def V_total_hat(R):
    """
    R: distance in Angstrom
    V: pair potential in Kelvin
    """
    V = V_BO_hat(R) + V_AD_hat(R) + V_REL_hat(R) + V_QED_hat(R)
    return V 

def dV_total_hat(R):
    """
    R: distance in Angstrom
    dV: derivative with respect to distance in Kelvin
    """        
    dV = derivative(V_total_hat,R,dx=1e-6)
    return dV

def U_tail(rc,plot=False):
    r = np.linspace(rc,1000,100000)
    dr = r[1]-r[0]
    U = V_total_hat(r)
    Uintegrand = U*r**2*dr
    Uint = Uintegrand.sum()
    
    if plot:
    
        plt.plot(r,U*r**2*dr,'k')
        plt.xlabel('Distance (A)')
        plt.ylabel(r'$U \times r^2 dr (K A^3)$')
        plt.show()
        
        Urun_sum = np.zeros(len(r))
        
        for i,U in enumerate(Uintegrand):
            
#            run_sum = (U[:i+1]*r[:i+1]**2*dr).sum()
            Urun_sum[i:] += U
            
        plt.plot(r,Urun_sum,'k')
        plt.xlabel('Cut-off Distance (A)')
        plt.ylabel(r'Tail Correction, U')
        plt.show()
        
        plt.plot(r,Urun_sum/Uint*100.,'k')
        plt.xlabel('Cut-off Distance (A)')
        plt.ylabel(r'Percent Tail Correction, U')
        plt.show()
        
    return Uint

def P_tail(rc,plot=False):
    r = np.linspace(rc,1000,100000)
    dr = r[1]-r[0]
    dU = dV_total_hat(r)
    Pintegrand = dU*r**3*dr
    Pint = Pintegrand.sum()
    
    if plot:
    
        plt.plot(r,dU*r**3*dr,'k')
        plt.xlabel('Distance (A)')
        plt.ylabel(r'$dU/dr \times r^3 dr (K A^3)$')
        plt.show()
        
        Prun_sum = np.zeros(len(r))
        
        for i,P in enumerate(Pintegrand):
            
#            run_sum = (U[:i+1]*r[:i+1]**2*dr).sum()
            Prun_sum[i:] += P
            
        plt.plot(r,Prun_sum,'k')
        plt.xlabel('Cut-off Distance (A)')
        plt.ylabel(r'Tail Correction, P')
        plt.show()
        
        plt.plot(r,Prun_sum/Pint*100.,'k')
        plt.xlabel('Cut-off Distance (A)')
        plt.ylabel(r'Percent Tail Correction, P')
        plt.show()
        
    return Pint

RV_BO_lit = np.loadtxt('RV_BO_lit.txt')
R_BO_lit = RV_BO_lit[:,0]*bohr_to_A
V_BO_lit = RV_BO_lit[:,1]

V_BO = V_BO_hat(R_BO_lit)

RV_AD_lit = np.loadtxt('RV_AD_lit.txt')
R_AD_lit = RV_AD_lit[:,0]*bohr_to_A
V_AD_lit = RV_AD_lit[:,1]

V_AD = V_AD_hat(R_AD_lit)

RV_REL_lit = np.loadtxt('RV_REL_lit.txt')
R_REL_lit = RV_REL_lit[:,0]*bohr_to_A
V_REL_lit = RV_REL_lit[:,1]

V_REL = V_REL_hat(R_REL_lit)

RV_QED_lit = np.loadtxt('RV_QED_lit.txt')
R_QED_lit = RV_QED_lit[:,0]*bohr_to_A
V_QED_lit = RV_QED_lit[:,1]

V_QED = V_QED_hat(R_QED_lit)

RV_total_lit = np.loadtxt('RV_total_lit.txt')
R_total_lit = RV_total_lit[:,0]*bohr_to_A
V_total_lit = RV_total_lit[:,1]

V_total = V_total_hat(R_total_lit)

R_plot = np.linspace(0.1,50,30000)*bohr_to_A
V_total_plot = V_total_hat(R_plot)
dV_total_plot = dV_total_hat(R_plot)

plt.plot(R_BO_lit,V_BO_lit,'bo',label='Literate')
plt.plot(R_BO_lit,V_BO,'b--',label='This Code')
plt.plot(R_AD_lit,V_AD_lit,'rx',label='Literate')
plt.plot(R_AD_lit,V_AD,'r--',label='This Code')
plt.plot(R_REL_lit,V_REL_lit,'cs',label='Literate')
plt.plot(R_REL_lit,V_REL,'c--',label='This Code')
plt.plot(R_QED_lit,V_QED_lit,'gv',label='Literate')
plt.plot(R_QED_lit,V_QED,'g--',label='This Code')
plt.plot(R_plot,V_total_plot,'k',label='This Code')
plt.xlabel('Distance (A)')
plt.ylabel('Pair Potential (K)')
plt.show()

plt.plot(R_BO_lit,V_BO_lit,'bo',label='Literate')
plt.plot(R_BO_lit,V_BO,'b--',label='This Code')
plt.plot(R_AD_lit,V_AD_lit,'rx',label='Literate')
plt.plot(R_AD_lit,V_AD,'r--',label='This Code')
plt.plot(R_REL_lit,V_REL_lit,'cs',label='Literate')
plt.plot(R_REL_lit,V_REL,'c--',label='This Code')
plt.plot(R_QED_lit,V_QED_lit,'gv',label='Literate')
plt.plot(R_QED_lit,V_QED,'g--',label='This Code')
plt.plot(R_plot,V_total_plot,'k',label='This Code')
plt.xlabel('Distance (A)')
plt.ylabel('Pair Potential (K)')
#plt.xlim([4*bohr_to_A,9*bohr_to_A])
plt.ylim([-20,100])
plt.show()

plt.plot(R_BO_lit,(V_BO - V_BO_lit)/V_BO_lit*100.,'bo')
plt.plot(R_AD_lit,(V_AD - V_AD_lit)/V_AD_lit*100.,'rx')
plt.plot(R_REL_lit,(V_REL - V_REL_lit)/V_REL_lit*100.,'cs')
plt.plot(R_QED_lit,(V_QED - V_QED_lit)/V_QED_lit*100.,'gv')
plt.plot(R_total_lit,(V_total - V_total_lit)/V_total_lit*100.,'ks')
plt.xlabel('Distance (A)')
plt.ylabel('Percent Deviation')
plt.show()

plt.plot(V_BO_lit,(V_BO - V_BO_lit)/V_BO_lit*100.,'bo')
plt.plot(V_AD_lit,(V_AD - V_AD_lit)/V_AD_lit*100.,'rx')
plt.plot(V_REL_lit,(V_REL - V_REL_lit)/V_REL_lit*100.,'cs')
plt.plot(V_QED_lit,(V_QED - V_QED_lit)/V_QED_lit*100.,'gv')
plt.plot(V_total_lit,(V_total - V_total_lit)/V_total_lit*100.,'ks')
plt.xlabel('Pair Potential (K)')
plt.ylabel('Percent Deviation')
plt.show()

plt.plot(V_BO_lit,V_BO - V_BO_lit,'bo')
plt.plot(V_AD_lit,V_AD - V_AD_lit,'rx')
plt.plot(V_REL_lit,V_REL - V_REL_lit,'cs')
plt.plot(V_QED_lit,V_QED - V_QED_lit,'gv')
plt.plot(V_total_lit,V_total - V_total_lit,'ks')
plt.xlabel('Pair Potential (K)')
plt.ylabel('Deviation (K)')
plt.show()

plt.plot(np.log10(np.abs(V_BO_lit)),(V_BO - V_BO_lit)/V_BO_lit*100.,'bo')
plt.plot(np.log10(np.abs(V_AD_lit)),(V_AD - V_AD_lit)/V_AD_lit*100.,'rx')
plt.plot(np.log10(np.abs(V_REL_lit)),(V_REL - V_REL_lit)/V_REL_lit*100.,'cs')
plt.plot(np.log10(np.abs(V_QED_lit)),(V_QED - V_QED_lit)/V_QED_lit*100.,'gv')
plt.plot(np.log10(np.abs(V_total_lit)),(V_total - V_total_lit)/V_total_lit*100.,'ks')
plt.xlabel('Log10 Absolute Value Pair Potential (K)')
plt.ylabel('Percent Deviation')
plt.show()

plt.plot(R_total_lit,V_total_lit,'ks',mfc='None',label='Literate')
plt.plot(R_plot,V_total_plot,'k',label='This Code')
plt.xlabel('Distance (A)')
plt.ylabel('Pair Potential (K)')
#plt.xlim([4*bohr_to_A,12*bohr_to_A])
plt.ylim([-20,100])
plt.show()

plt.semilogy(R_total_lit,V_total_lit,'ks',mfc='None',label='Literate')
plt.semilogy(R_plot,V_total_plot,'k',label='This Code')
plt.xlabel('Distance (A)')
plt.ylabel('Pair Potential (K)')
plt.xlim([0,4*bohr_to_A])
plt.ylim([100,10000000])
plt.show()

plt.plot(R_total_lit,V_total_lit,'ks',mfc='None',label='Literate')
plt.plot(R_plot,V_total_plot,'k',label='This Code')
plt.xlabel('Distance (A)')
plt.ylabel('Pair Potential (K)')
plt.xlim([10*bohr_to_A,50*bohr_to_A])
plt.ylim([-0.5,0.005])
plt.show()

plt.plot(V_total_lit,(V_total - V_total_lit)/V_total_lit*100.,'ks')
plt.xlabel('Pair Potential (K)')
plt.ylabel('Percent Deviation')
plt.show()

plt.plot(V_total_lit,V_total - V_total_lit,'ks')
plt.xlabel('Pair Potential (K)')
plt.ylabel('Deviation (K)')
plt.show()

plt.plot(R_plot,-dV_total_plot,'k--')
plt.xlim([0.5,14])
plt.ylim([-20,100])
plt.xlabel('Distance (A)')
plt.ylabel('Force (K/A)')
plt.show()

#R_towhee = np.linspace(0.05,14,3000)
R_low = np.linspace(0.05,1.5,5)
R_mid = np.linspace(R_low[-1]+R_low[-1]-R_low[-2],10,2495)
R_high = np.linspace(R_mid[-1]+R_mid[-1]-R_mid[-2],14,500)
R_towhee = np.concatenate([R_low,R_mid,R_high])
V_towhee = V_total_hat(R_towhee)
V_towhee[0:4] = V_towhee.max() 

print(U_tail(14,True))
print(P_tail(14,True))

### Fitting the pair potential to models that are in Towhee:
    
rmin = 0.265 #[nm]
rmax = 0.6  #[nm]

def double_exponential(params,r):
    
    if params[1] > -1e5 and params[3] > -1e5:
        U = params[0]*np.exp(-params[1]*r) - params[2]*np.exp(-params[3]*r)
    else:
        U = 0.
    return U    
    
guess = np.array([-1418.200759,15.37324216,-84995699.95,56.97745683])

rfit = R_towhee/10. #[nm]
Vfit = V_towhee[rfit>rmin]
rfit = rfit[rfit>rmin]
Vfit = Vfit[rfit<rmax]
rfit = rfit[rfit<rmax]

SSE = lambda params: np.sum((Vfit - double_exponential(params,rfit))**2)

SSE(guess) 

opt = minimize(SSE,guess)  
#opt = basinhopping(SSE,guess,niter=1000) #This did not give me a better result when I gave it a very different guess

plt.plot(R_towhee/10.,V_towhee,'k-')
plt.plot(R_towhee/10.,double_exponential(opt.x,R_towhee/10.),'r--')
plt.xlim([0.2,0.6])
plt.ylim([-12,12])
plt.show()

plt.plot(R_towhee/10.,V_towhee,'k-')
plt.plot(R_towhee/10.,double_exponential(opt.x,R_towhee/10.),'r--')
plt.xlim([0.6,1.4])
plt.ylim([-0.3,0])
plt.show()

plt.plot(R_towhee/10.,V_towhee,'k-')
plt.plot(R_towhee/10.,double_exponential(opt.x,R_towhee/10.),'r--')
plt.xlim([0.1,0.2])
plt.ylim([10,10000])
plt.show()
