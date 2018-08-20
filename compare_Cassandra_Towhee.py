# -*- coding: utf-8 -*-
"""
Plot Helium results
"""

from __future__ import division
import numpy as np 
import os, sys, argparse, shutil
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

#sys.path.insert(0,'H:/Python/Lennard-Jones/')
#
#from LJ_fluid_correlations import *

font = {'size' : '24'}
plt.rc('font',**font)

root_path = 'H:/Helium_ab_initio/'

Nmol = 2800
rcut = 14 #[Ang]

packages = ['Cassandra','Towhee','Cassandra_LJ','Towhee_LJ_native','Towhee_LJ_tabulated']   

color_scheme = {'Cassandra':'b','Towhee':'r','Cassandra_LJ':'g','Towhee_LJ_native':'c','Towhee_LJ_tabulated':'m'}
shape_scheme = {'Cassandra':'o','Towhee':'s','Cassandra_LJ':'^','Towhee_LJ_native':'v','Towhee_LJ_tabulated':'d'}

fig = plt.figure(figsize=(12,12)) 

for pack in packages:

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
    Temp_all = VLCC[:,0]
    rho_liq = VLCC[:,1]
    rho_vap = VLCC[:,2]
    
    rho_rect = (rho_liq+rho_vap)/2.

    plt.plot(rho_liq,Temp_all,color_scheme[pack]+shape_scheme[pack],markersize=8,mfc='None',label=pack)
    plt.plot(rho_vap,Temp_all,color_scheme[pack]+shape_scheme[pack],markersize=8,mfc='None')
    plt.plot(rho_rect,Temp_all,color_scheme[pack]+shape_scheme[pack],markersize=8,mfc='None')
    
### Copied from LJ_fluid_correlations
# Conversion constants

k_B = 1.38065e-23 #[J/K]
N_A = 6.02214e23 #[1/mol]
m3_to_nm3 = 1e27
gm_to_kg = 1./1000
J_to_kJ = 1./1000
J_per_m3_to_kPA = 1./1000

T_c_star= 1.31
rho_c_star= 0.314
rho_L_star_params= [0.477, 0.2124, -0.01151]
rho_v_star_params= [-0.477, 0.05333, 0.1261]
P_v_star_params = [1.2629, -4.8095, -0.15115]

rho_L_star_params = [rho_c_star, T_c_star, rho_L_star_params[0], rho_L_star_params[1], rho_L_star_params[2]]
rho_v_star_params = [rho_c_star, T_c_star, rho_v_star_params[0], rho_v_star_params[1], rho_v_star_params[2]]

def rho_L_star_hat(T_star, b = rho_L_star_params):
    tau = np.ones(len(T_star))*b[1] - T_star # T_c_star - T_star
    rho_L_star = b[0] + b[2]*tau**(1./3) + b[3]*tau + b[4]*tau**(3./2)
    return rho_L_star

def rho_L_hat_LJ(T,eps,sig,M_w):
    T_star = T/(np.ones(len(T))*eps)
    rho_L_star = rho_L_star_hat(T_star)
    rho_L = rho_L_star *  M_w  / sig**3 / N_A * m3_to_nm3 * gm_to_kg #[kg/m3]
    return rho_L

def rho_v_star_hat(T_star, b = rho_v_star_params):
    tau = np.ones(len(T_star))*b[1] - T_star # T_c_star - T_star
    rho_v_star = b[0] + b[2]*tau**(1./3) + b[3]*tau + b[4]*tau**(3./2)
    return rho_v_star

def rho_v_hat_LJ(T,eps,sig,M_w):
    T_star = T/(np.ones(len(T))*eps)
    rho_v_star = rho_v_star_hat(T_star)
    rho_v = rho_v_star *  M_w  / sig**3 / N_A * m3_to_nm3 * gm_to_kg #[kg/m3]
    return rho_v

def P_v_star_hat(T_star, b = P_v_star_params):
    P_v_star = np.exp(b[0]*T_star + b[1]/T_star + b[2]/(T_star**4))
    return P_v_star

def P_v_hat_LJ(T,eps,sig):
    T_star = T/(np.ones(len(T))*eps)
    P_v_star = P_v_star_hat(T_star)
    P_v = P_v_star *  eps  / sig**3 * k_B * m3_to_nm3 * J_per_m3_to_kPA #[kPa]
    return P_v

### Properties for Helium
sig_LJ = 0.264 #[nm]
eps_LJ = 11. #[K]
M_w = 4.0026 #[gm/mol]

T_plot = np.linspace(np.min(Temp_all),T_c_star*eps_LJ,500)
rho_L_plot = rho_L_hat_LJ(T_plot,eps_LJ,sig_LJ,M_w)
rho_v_plot = rho_v_hat_LJ(T_plot,eps_LJ,sig_LJ,M_w)
P_v_plot = P_v_hat_LJ(T_plot,eps_LJ,sig_LJ)/100. #[bar]

plt.plot(rho_L_plot,T_plot,'g--',label='Lennard-Jones')
plt.plot(rho_v_plot,T_plot,'g--')
#plt.plot([311.,250.],[7.,11.],'go',markersize=10,mfc='None',label='Cassandra, Lennard-Jones')

plt.ylabel('Temperature (K)')
plt.xlabel('Density (kg/m3)')
plt.xlim([-10,None])
plt.ylim([None,12])
plt.yticks([7,8,9,10,11,12])
plt.tight_layout()
plt.legend()
fig.savefig(root_path+'VLCC_Cassandra_Towhee.pdf')
plt.close()  
