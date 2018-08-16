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

root_path = 'H:/Helium_ab_initio/Towhee/Results/'

Nmol = 2800
rcut = 14 #[Ang]

system_path = root_path+str(Nmol)+'_'+str(rcut)+'/'

Temps = np.array([7.,8.,9.,10.,11.])                                     
                                     
NTemps = len(Temps)
NReps=8
NBoxes=2
NSims = NTemps*NReps
MCSteps = 100000

nStage=2

rho_liq = np.zeros(NSims)
rho_vap = np.zeros(NSims)
press_vap = np.zeros(NSims)
Temp_all = np.zeros(NSims)

counter = 0

for iTemp, Temp_i in enumerate(Temps):
    
    Temp_path = system_path+str(int(Temp_i))+'_'
    
    for iRep in range(NReps):                            
                                                                
        rep_path = Temp_path+str(iRep)+'/'
        
#        if os.path.exists(rep_path+'box_01_step_00000000100000.pdb'):
#            
#            print('file exists')
#            
#            for iBox in np.arange(1,NBoxes+1):
#                
#                towhee_vlcc = np.loadtxt(rep_path+'towhee_vlcc')
#                
#                rho_liq[counter] = towhee_vlcc[3]
#                rho_vap[counter] = towhee_vlcc[1]
#                press_vap[counter] = towhee_vlcc[7]
#            
#            Temp_all[counter] = Temp_i
#            counter += 1
#            
#        else:
#            
#            Output_2 = np.loadtxt(rep_path+'Output_2',skiprows=1152)
#            print(Output_2)

        rho_box_1 = np.loadtxt(rep_path+'box_1_averages_'+str(nStage))
        rho_box_2 = np.loadtxt(rep_path+'box_2_averages_'+str(nStage))    
        
        rho_liq_blocks = np.max([rho_box_1,rho_box_2],axis=0)
        rho_vap_blocks = np.min([rho_box_1,rho_box_2],axis=0)
        
        rho_liq[counter] = np.mean(rho_liq_blocks)
        rho_vap[counter] = np.mean(rho_vap_blocks)           

        if np.mean(rho_box_1) != np.mean(rho_liq_blocks):

            print('WARNING: Vapor and liquid phase have swapped boxes')
        
        Temp_all[counter] = Temp_i
        counter += 1
            
rho_liq = rho_liq[Temp_all>0]
rho_vap = rho_vap[Temp_all>0]
press_vap = press_vap[Temp_all>0]
Temp_all = Temp_all[Temp_all>0]
       
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

if False:

    fig = plt.figure(figsize=(12,12))
    plt.plot(1000./Temp_all,np.log10(press_vap),'bo',markersize=5,mfc='None')
    plt.xlabel(r'1000/Temperature (K)')
    plt.ylabel(r'log$_{10}$(Pressure/bar)')
    plt.tight_layout()
    fig.savefig(system_path+'plots/ClausiusClapeyron.pdf')
    plt.close()  
    
f = open(system_path+'VLCC.txt','w')

f.write(r'Tsat (K)'+'\t'+r'rho_liq (kg/m3)'+'\t'+r'rho_vap (kg/m3)')

for Tsat_i, rho_liq_i, rho_vap_i in zip(Temp_all,rho_liq,rho_vap):
    
    f.write('\n'+str(Tsat_i)+'\t'+str(rho_liq_i*1000.)+'\t'+str(rho_vap_i*1000.))
    
f.close()