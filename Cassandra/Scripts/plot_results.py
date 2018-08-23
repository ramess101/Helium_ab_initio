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

root_path = '/home/ram9/Cassandra_Helium/'

Nmol = 1400
rcut = 14 #[Ang]

system_path = root_path+str(Nmol)+'_'+str(rcut)+'/'

Temps = np.array([7.,8.,9.,10.,11.])                                     
                                     
NTemps = len(Temps)
NReps=3
NBoxes=2
NSims = NTemps*NReps

#fig_VLCC = plt.figure(figsize=(6,6))
#fig_Antoine = plt.figure(figsize=(6,6))
#
#plt.close()

rho_liq = np.zeros(NSims)
rho_vap = np.zeros(NSims)
press_vap = np.zeros(NSims)
Temp_all = np.zeros(NSims)

counter = 0

for iTemp, Temp_i in enumerate(Temps):
    
    Temp_path = system_path+'temp_'+str(iTemp)+'/'
    
    for iRep in range(NReps):                            
                                                                
        rep_path = Temp_path+'rep_'+str(iRep)+'/'
        
        for iBox in np.arange(1,NBoxes+1):
            
            prp_box = np.loadtxt(rep_path+'helium.out.box'+str(iBox)+'.prp',skiprows=3)
            
            MC_step = prp_box[:,0]
            
            Nstep_equil = 2000000
            Nstep_all = MC_step[-1]

            rho_box_all = prp_box[:,2] #[kg/m3]
            rho_box_prod = rho_box_all[MC_step>Nstep_equil]
            rho_box_avg = np.mean(rho_box_prod)            

            press_box_all = prp_box[:,4] # [bar]
            press_box_prod = press_box_all[MC_step>Nstep_equil]
            press_box_avg = np.mean(press_box_prod)            

#            fig = plt.figure(figsize=(6,6))
 #           plt.plot(MC_step,rho_box_all,'k-')
  #          plt.plot([Nstep_equil,Nstep_all],[rho_box_avg,rho_box_avg],'r--')
   #         plt.xlabel('MC step')
    #        plt.ylabel('Density (kg/m3)')
     #       plt.tight_layout()
      #      fig.savefig(system_path+'plots/temp_'+str(iTemp)+'_rep_'+str(iRep)+'_box_'+str(iBox)+'_rho.pdf')
       #     plt.close()  
            
#            fig = plt.figure(figsize=(6,6))
 #           plt.plot(MC_step,press_box_all,'k-')
  #          plt.plot([Nstep_equil,Nstep_all],[press_box_avg,press_box_avg],'r--')
   #         plt.xlabel('MC step')
    #        plt.ylabel('Pressure (bar)')
     #       plt.tight_layout()
      #      fig.savefig(system_path+'plots/temp_'+str(iTemp)+'_rep_'+str(iRep)+'_box_'+str(iBox)+'_press.pdf')
       #     plt.close()

            if iBox == 1:

                rho_liq[counter] = rho_box_avg
                rho_box_1 = rho_box_prod                

            if iBox == 2:

                rho_vap[counter] = rho_box_avg
                press_vap[counter] = press_box_avg
            
                for rho_box_1_step, rho_box_2_step in zip(rho_box_1,rho_box_prod):

                    if rho_box_2_step > rho_box_1_step:

                        print('WARNING: Vapor and liquid phase have swapped boxes')
        
        Temp_all[counter] = Temp_i
        counter += 1
       

#print(rho_liq)

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

f = open(system_path+'VLCC.txt','w')

f.write(r'Tsat (K)'+'\t'+r'rho_liq (kg/m3)'+'\t'+r'rho_vap (kg/m3)')

for Tsat_i, rho_liq_i, rho_vap_i in zip(Temp_all,rho_liq,rho_vap):
    
    f.write('\n'+str(Tsat_i)+'\t'+str(rho_liq_i)+'\t'+str(rho_vap_i))
    
f.close()

f = open(system_path+'Psat.txt','w')

f.write(r'Tsat (K)'+'\t'+r'Psat (bar)')

for Tsat_i, press_vap_i in zip(Temp_all,press_vap):
    
    f.write('\n'+str(Tsat_i)+'\t'+str(press_vap_i))
    
f.close()