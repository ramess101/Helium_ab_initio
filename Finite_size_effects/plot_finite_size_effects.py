# -*- coding: utf-8 -*-
"""
Plot Helium results
"""

from __future__ import division
import numpy as np 
import os, sys, argparse, shutil
import matplotlib.pyplot as plt

font = {'size' : '24'}
plt.rc('font',**font)

def compile_VLCC(root_path,Temps,nreps):

    f0 = open(root_path+'VLCC.txt','w')
    f0.write('Tsat (K)'+'\t'+'rhol (kg/m3)'+'\t'+'rhov (kg/m3)'+'\n')
    
    for Tsat in Temps:
        for irep in range(0,nreps):
            box1_all = np.loadtxt(root_path+'/'+str(Tsat)+'_'+str(irep)+'/box_1_averages_2')
            box2_all = np.loadtxt(root_path+'/'+str(Tsat)+'_'+str(irep)+'/box_2_averages_2')    
            box1_avg = np.mean(box1_all)*1000
            box2_avg = np.mean(box2_all)*1000
            
            f0.write(str(Tsat)+'\t'+str(box1_avg)+'\t'+str(box2_avg)+'\n')
            
    f0.close()

root_path = 'H:/Helium_ab_initio/'

system_list = ['Towhee_2800_18','Cassandra_2800_14','Cassandra_1400_10','Cassandra_800_10']

Nmol_dic = {'Towhee_2800_14':2800,'Towhee_2800_18':2800,'Cassandra_1400_10':1400,'Cassandra_2800_14':2800,'Cassandra_800_10':800}
rcut_dic = {'Towhee_2800_14':1.4,'Towhee_2800_18':1.8,'Cassandra_1400_10':1.0,'Cassandra_2800_14':1.4,'Cassandra_800_10':1.0}

color_dic = {'Towhee_2800_14':'m','Towhee_2800_18':'g','Cassandra_1400_10':'r','Cassandra_2800_14':'b','Cassandra_800_10':'c'}
shape_dic = {'Towhee_2800_14':'^','Towhee_2800_18':'s','Cassandra_1400_10':'d','Cassandra_2800_14':'o','Cassandra_800_10':'v'}

tail_path_dic = {'Towhee_2800_14':'Towhee/Results/2800_14/','Towhee_2800_18':'Towhee/Results/2800_18/','Cassandra_1400_10':'Cassandra/Results/1400_10_md1/','Cassandra_2800_14':'Cassandra/Results/Helium_2800_14/','Cassandra_800_10':'Cassandra/Results/3Body/800_10_all/'}

label_dic = {'Towhee_2800_14':r'$N = 2800$, $r_{\rm cut} = 1.4$ nm','Towhee_2800_18':r'$N = 2800$, $r_{\rm cut} = 1.8$ nm','Cassandra_1400_10':r'$N = 1400$, $r_{\rm cut} = 1.0$ nm','Cassandra_2800_14':r'$N = 2800$, $r_{\rm cut} = 1.4$ nm','Cassandra_800_10':r'$N = 800$, $r_{\rm cut} = 1.0$ nm'}

##    #### Limit the range of data included in the fit
TsatLow = 7
TsatHigh = 11.1
Tplot_low = 9

Mw_He = 4.0026 #[gm/mol]

fig0, ax0 = plt.subplots(nrows=1,ncols=1,figsize=(10,10)) 

for system in system_list:
    
    Nmol = Nmol_dic[system]
    rcut = rcut_dic[system]
    color = color_dic[system]
    shape = shape_dic[system]
    system_path = root_path+tail_path_dic[system] 
    label = label_dic[system]                             

    try:
        VLCC = np.loadtxt(system_path+'VLCC.txt',skiprows=1)
    except:
        compile_VLCC(system_path,Temps=[7,8,9,10,11],nreps=8)
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
    
ax0.set_ylim([8.9,11.1])
ax0.set_yticks([9.0,9.5,10.0,10.5,11.0])
#ax0.set_ylim([5.5,11.1])
#ax0.set_yticks([7,8,9,10,11])
ax0.set_xlabel(r'$\rho$ (kmol/m$^3$)',fontsize=28)    
ax0.set_ylabel('$T$ (K)',fontsize=28)
#ax0.legend()

lgd = ax0.legend(loc='lower center', bbox_to_anchor=(0.43,-0.27),
          ncol=2,numpoints=1,handlelength=2,handletextpad=0.2,columnspacing=0.5,frameon=True,borderaxespad=0)

plt.tight_layout()

fig0.savefig('GEMC_finite_size_effects.pdf',bbox_extra_artists=(lgd,),bbox_inches='tight')