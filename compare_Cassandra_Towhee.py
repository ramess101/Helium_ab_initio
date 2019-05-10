# -*- coding: utf-8 -*-
"""
Plot Helium results
"""

from __future__ import division
import numpy as np 
import os, sys, argparse, shutil
#import matplotlib
#matplotlib.use('Agg')
import matplotlib.pyplot as plt

#sys.path.insert(0,'H:/Python/Lennard-Jones/')
#
#from LJ_fluid_correlations import *

font = {'size' : '24'}
plt.rc('font',**font)

root_path = 'C:/Users/rmesserl/Documents/NIST_projects/Helium_ab_initio/'

Nmol = 2800
rcut = 14 #[Ang]

packages = ['Cassandra','Towhee'] #,'Cassandra_LJ','Towhee_LJ_native','Towhee_LJ_tabulated']   

color_dic = {'Cassandra':'b','Towhee':'r','Cassandra_LJ':'g','Towhee_LJ_native':'c','Towhee_LJ_tabulated':'m'}
shape_dic = {'Cassandra':'o','Towhee':'s','Cassandra_LJ':'^','Towhee_LJ_native':'v','Towhee_LJ_tabulated':'d'}
label_dic = {'Cassandra':r'Cassandra','Towhee':r'Towhee'}

#fig, ax = plt.subplots(nrows=1,ncols=1,figsize=(10,10)) 
#ax.set_xscale('log')

##    #### Limit the range of data included in the fit
TsatLow = 7
TsatHigh = 15
Tplot_low = 7

Mw_He = 4.0026 #[gm/mol]

fig0, ax0 = plt.subplots(nrows=1,ncols=1,figsize=(10,10)) 

for pack in packages:
    
    color = color_dic[pack]
    shape = shape_dic[pack]
    label = label_dic[pack]

    if pack == 'Cassandra':
    
        system_path = root_path+pack+'/Results/Helium_'+str(Nmol)+'_'+str(rcut)+'/'                                    
    
    elif pack == 'Towhee':
    
        system_path = root_path+pack+'/Results/'+str(Nmol)+'_'+str(rcut)+'/'                                    

    elif pack == 'Cassandra_LJ':
    
        system_path = root_path+'Cassandra/Results/LJ_1400_14/'                                    

    elif pack == 'Towhee_LJ':
    
        system_path = root_path+'Towhee/Results/LJ_truncated_native_2800_18/'
        
    elif pack == 'Towhee_LJ_tabulated':
    
        system_path = root_path+'Towhee/Results/LJ_truncated_2800_18/'

    VLCC = np.loadtxt(system_path+'VLCC.txt',skiprows=1)
    Tsat = VLCC[:,0]
    rhol = VLCC[:,1]
    rhov = VLCC[:,2]
               
    ### Trim low and/or high data
        
    rhol = rhol[Tsat>=TsatLow]
    rhov = rhov[Tsat>=TsatLow]
    Tsat = Tsat[Tsat>=TsatLow]
    
    rhol = rhol[Tsat<=TsatHigh]
    rhov = rhov[Tsat<=TsatHigh]
    Tsat = Tsat[Tsat<=TsatHigh]
    
    rhor = (rhol+rhov)/2.
    
    Tsat_avg = np.unique(Tsat)
    
    rhol_avg = np.zeros(len(Tsat_avg))
    rhov_avg = np.zeros(len(Tsat_avg))
    rhor_avg = np.zeros(len(Tsat_avg))
    rhol_95 = np.zeros(len(Tsat_avg))
    rhov_95 = np.zeros(len(Tsat_avg))
    rhor_95 = np.zeros(len(Tsat_avg))
    
    for iTsat, Tsat_i in enumerate(Tsat_avg):
        rhol_avg[iTsat] = np.mean(rhol[Tsat == Tsat_i])
        rhov_avg[iTsat] = np.mean(rhov[Tsat == Tsat_i])
        rhor_avg[iTsat] = np.mean(rhor[Tsat == Tsat_i])

        rhol_95[iTsat] = 3.18 * np.std(rhol[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        rhov_95[iTsat] = 3.18 * np.std(rhov[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        rhor_95[iTsat] = 3.18 * np.std(rhor[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        
    ax0.errorbar(rhol_avg[Tsat_avg >= Tplot_low]/Mw_He,Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhol_95[Tsat_avg >= Tplot_low]/Mw_He,fmt=color+shape,mfc='None',markersize=10,label=label,capsize=6)
    ax0.errorbar(rhov_avg[Tsat_avg >= Tplot_low]/Mw_He,Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhov_95[Tsat_avg >= Tplot_low]/Mw_He,fmt=color+shape,mfc='None',markersize=10,capsize=6)
    ax0.errorbar(rhor_avg[Tsat_avg >= Tplot_low]/Mw_He,Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhor_95[Tsat_avg >= Tplot_low]/Mw_He,fmt=color+shape,mfc='None',markersize=10,capsize=6)  
    
#ax0.set_ylim([6.9,11.1])
#ax0.set_yticks([7,8,9,10,11])
ax0.set_ylim([8.9,11.1])
ax0.set_yticks([9,10,11])
ax0.set_xlabel(r'$\rho$ (mol$\cdot$L$^{-1}$)',fontsize=28)    
ax0.set_ylabel('$T$ (K)',fontsize=28)

#ax0.legend()
plt.tight_layout()
#
#fig0.savefig(root_path+'VLCC_Cassandra_Towhee.pdf')

lgd = ax0.legend(loc='lower center', bbox_to_anchor=(0.5,-0.22),
          ncol=2,numpoints=1,handlelength=2,handletextpad=0.2,columnspacing=0.5,frameon=True,borderaxespad=0)

plt.tight_layout()

fig0.savefig(root_path+'VLCC_Cassandra_Towhee.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')

#    plt.plot(rho_liq,Temp_all,color_scheme[pack]+shape_scheme[pack],markersize=8,mfc='None',label=pack)
#    plt.plot(rho_vap,Temp_all,color_scheme[pack]+shape_scheme[pack],markersize=8,mfc='None')
#    plt.plot(rho_rect,Temp_all,color_scheme[pack]+shape_scheme[pack],markersize=8,mfc='None')
    
#    if pack == 'Cassandra':
#        
#        rho_liq_Cass = rho_liq.copy()
#        Temp_Cass = Temp_all.copy()
#        
#    elif pack == 'Towhee':
#        
#        rho_liq_Tow = rho_liq.copy()
#        Temp_Tow = Temp_all.copy()
#        
#for Temp in np.unique(Temp_Tow):
#    print(np.mean(rho_liq_Cass[Temp_Cass==Temp]/rho_liq_Tow[Temp_Tow==Temp]))
    
### Copied from LJ_fluid_correlations
# Conversion constants
#
#k_B = 1.38065e-23 #[J/K]
#N_A = 6.02214e23 #[1/mol]
#m3_to_nm3 = 1e27
#gm_to_kg = 1./1000
#J_to_kJ = 1./1000
#J_per_m3_to_kPA = 1./1000
#
#T_c_star= 1.31
#rho_c_star= 0.314
#rho_L_star_params= [0.477, 0.2124, -0.01151]
#rho_v_star_params= [-0.477, 0.05333, 0.1261]
#P_v_star_params = [1.2629, -4.8095, -0.15115]
#
#rho_L_star_params = [rho_c_star, T_c_star, rho_L_star_params[0], rho_L_star_params[1], rho_L_star_params[2]]
#rho_v_star_params = [rho_c_star, T_c_star, rho_v_star_params[0], rho_v_star_params[1], rho_v_star_params[2]]
#
#def rho_L_star_hat(T_star, b = rho_L_star_params):
#    tau = np.ones(len(T_star))*b[1] - T_star # T_c_star - T_star
#    rho_L_star = b[0] + b[2]*tau**(1./3) + b[3]*tau + b[4]*tau**(3./2)
#    return rho_L_star
#
#def rho_L_hat_LJ(T,eps,sig,M_w):
#    T_star = T/(np.ones(len(T))*eps)
#    rho_L_star = rho_L_star_hat(T_star)
#    rho_L = rho_L_star *  M_w  / sig**3 / N_A * m3_to_nm3 * gm_to_kg #[kg/m3]
#    return rho_L
#
#def rho_v_star_hat(T_star, b = rho_v_star_params):
#    tau = np.ones(len(T_star))*b[1] - T_star # T_c_star - T_star
#    rho_v_star = b[0] + b[2]*tau**(1./3) + b[3]*tau + b[4]*tau**(3./2)
#    return rho_v_star
#
#def rho_v_hat_LJ(T,eps,sig,M_w):
#    T_star = T/(np.ones(len(T))*eps)
#    rho_v_star = rho_v_star_hat(T_star)
#    rho_v = rho_v_star *  M_w  / sig**3 / N_A * m3_to_nm3 * gm_to_kg #[kg/m3]
#    return rho_v
#
#def P_v_star_hat(T_star, b = P_v_star_params):
#    P_v_star = np.exp(b[0]*T_star + b[1]/T_star + b[2]/(T_star**4))
#    return P_v_star
#
#def P_v_hat_LJ(T,eps,sig):
#    T_star = T/(np.ones(len(T))*eps)
#    P_v_star = P_v_star_hat(T_star)
#    P_v = P_v_star *  eps  / sig**3 * k_B * m3_to_nm3 * J_per_m3_to_kPA #[kPa]
#    return P_v
#
#### Properties for Helium
#sig_LJ = 0.264 #[nm]
#eps_LJ = 11. #[K]
#M_w = 4.0026 #[gm/mol]
#
#T_plot = np.linspace(np.min(Temp_all),T_c_star*eps_LJ,500)
#rho_L_plot = rho_L_hat_LJ(T_plot,eps_LJ,sig_LJ,M_w)
#rho_v_plot = rho_v_hat_LJ(T_plot,eps_LJ,sig_LJ,M_w)
#P_v_plot = P_v_hat_LJ(T_plot,eps_LJ,sig_LJ)/100. #[bar]
#
##plt.plot(rho_L_plot,T_plot,'g--',label='Lennard-Jones')
##plt.plot(rho_v_plot,T_plot,'g--')
##plt.plot([311.,250.],[7.,11.],'go',markersize=10,mfc='None',label='Cassandra, Lennard-Jones')

#plt.ylabel('Temperature (K)')
#plt.xlabel('Density (kg/m3)')
#plt.xlim([-10,None])
#plt.ylim([None,12])
#plt.yticks([7,8,9,10,11,12])
#plt.tight_layout()
#plt.legend()
#fig.savefig(root_path+'VLCC_Cassandra_Towhee.pdf')
#plt.close()  
