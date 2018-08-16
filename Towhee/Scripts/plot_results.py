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

font = {'size' : '24'}
plt.rc('font',**font)

root_path = '/data/ram9/Helium/2Body_Total/'

Nmol = 2800
rcut = 14 #[Ang]

system_path = root_path+str(Nmol)+'_'+str(rcut)+'/'

Temps = np.array([7.,8.,9.,10.,11.])                                     
                                     
NTemps = len(Temps)
NReps=1
NBoxes=2
NSims = NTemps*NReps

rho_liq = np.zeros(NSims)
rho_vap = np.zeros(NSims)
press_vap = np.zeros(NSims)
Temp_all = np.zeros(NSims)

counter = 0

for iTemp, Temp_i in enumerate(Temps):
    
    Temp_path = system_path+str(Temp_i)+'_'
    
    for iRep in range(NReps):                            
                                                                
        rep_path = Temp_path+str(iRep)+'/'
        
        for iBox in np.arange(1,NBoxes+1):
            
            towhee_vlcc = np.loadtxt(rep_path+'towhee_vlcc')
            
            rho_liq[counter] = towhee_vlcc[3]
            rho_vap[counter] = towhee_vlcc[1]
            press_vap[counter] = towhee_vlcc[7]
        
        Temp_all[counter] = Temp_i
        counter += 1
       
rho_rect = (rho_liq + rho_vap)/2.

fig = plt.figure(figsize=(12,12))
plt.plot(rho_liq,Temp_all,'bo',markersize=5,mfc='None')
plt.plot(rho_vap,Temp_all,'go',markersize=5,mfc='None')
plt.plot(rho_rect,Temp_all,'co',markersize=5,mfc='None')
plt.ylabel('Temperature (K)')
plt.xlabel('Density (kg/m3)')
plt.tight_layout()
fig.savefig(system_path+'plots/VLCC.pdf')
plt.close()  

fig = plt.figure(figsize=(12,12))
plt.plot(1000./Temp_all,np.log10(press_vap),'bo',markersize=5,mfc='None')
plt.xlabel(r'1000/Temperature (K)')
plt.ylabel(r'log$_{10}$(Pressure/bar)')
plt.tight_layout()
fig.savefig(system_path+'plots/ClausiusClapeyron.pdf')
plt.close()  