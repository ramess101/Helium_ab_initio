# -*- coding: utf-8 -*-
"""
Creates publication plots for GEMC results
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy import stats
from scipy.optimize import minimize  
from VLE_model_fit import VLE

font = {'size' : '24'}
plt.rc('font',**font)

Mw_He = 4.0026 #[gm/mol]
Rg = 8.3144598e-5 #[m3 bar / K / mol]
         
Pc_Kofke = 0.928*10. #[bar]
rhoc_Kofke = 28.3 * Mw_He #[kg/m3]
Tc_Kofke = 13.00 #[K]

uPc_Kofke = 0.013*10. #[bar]
urhoc_Kofke = 0.4 * Mw_He #[kg/m3]
uTc_Kofke = 0.04 #[K]    

ulogPc_Kofke = (np.log10(Pc_Kofke+uPc_Kofke) - np.log10(Pc_Kofke-uPc_Kofke))/2.
uinvTc_Kofke = (10./(Tc_Kofke-uTc_Kofke)-10./(Tc_Kofke+uTc_Kofke))/2. 
               
Vc_Kofke = Mw_He/rhoc_Kofke/1000.
Zc_Kofke = Pc_Kofke * Vc_Kofke / Rg / Tc_Kofke

uZc_Kofke = Zc_Kofke * np.sqrt((uPc_Kofke/Pc_Kofke)**2. + (urhoc_Kofke/rhoc_Kofke)**2. + (uTc_Kofke/Tc_Kofke)**2.)
               
systems = ['1400_10_2Body','3Body_800_10']#'3Body_1400_10']#'3Body_800_10'] #,'2800_14_2Body']
path_dic = {'1400_10_2Body':'1400_10_md1/','2800_14_2Body':'Helium_2800_14/','3Body_800_10':'3Body/800_10_all/','3Body_1400_10':'3Body/1400_10_all/'}
label_dic = {'1400_10_2Body':'2-body','2800_14_2Body':'2-body','3Body_800_10':'(2+3)-body','3Body_1400_10':'(2+3)-body'}
color_dic = {'1400_10_2Body':'r','2800_14_2Body':'b','3Body_800_10':'b','3Body_1400_10':'b'}
shape_dic = {'1400_10_2Body':'o','2800_14_2Body':'s','3Body_800_10':'s','3Body_1400_10':'s'}
line_dic = {'1400_10_2Body':'-','2800_14_2Body':'--','3Body_800_10':'--','3Body_1400_10':'--'}

path_root = 'C:/Users/rmesserl/Documents/Helium_ab_initio/Cassandra/Results/'

##    #### Limit the range of data included in the fit
TsatLow = 9
TsatHigh = 15

fig0,ax0 = plt.subplots(ncols=1,nrows=2,figsize=[10,20])

ax0[0].set_xlabel(r'$\rho$ (mol$\cdot$L$^{-1}$)',fontsize=28)    
ax0[0].set_ylabel('$T$ (K)',fontsize=28)
ax0[1].set_xlabel(r'10/$T$ (1/K)',fontsize=28)
ax0[1].set_ylabel(r'$\log_{10}$($P$/MPa)',fontsize=28)

plt.tight_layout()

fig1,ax1 = plt.subplots(ncols=1,nrows=1,figsize=[10,10])

ax1.set_ylabel(r'$Z_{\rm vap}^{\rm sat}$')
ax1.set_xlabel(r'$T$ (K)')

plt.tight_layout()

#plt.tight_layout(pad=2)

for system in systems:
    
    path_tail = path_dic[system]
    system_label = label_dic[system]
    color = color_dic[system]
    shape = shape_dic[system]
    line = line_dic[system]
    
    #Load data
    path = path_root+path_tail
    VLCC = np.loadtxt(path+'VLCC.txt',skiprows=1)
    Psat = np.loadtxt(path+'Psat.txt',skiprows=1)[:,1]
    Tsat = VLCC[:,0]
    rhol = VLCC[:,1]
    rhov = VLCC[:,2]
    
    ### Trim low and/or high data
        
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

    rhor = (rhol+rhov)/2.
    logPsat = np.log10(Psat)
    Tsat_avg = np.unique(Tsat)
    invTsat_avg = 10./Tsat_avg
    
    Vv = Mw_He/rhov/1000. #[m3/mol]
    Zv = Psat*Vv/Rg/Tsat
    
    Psat_avg = np.zeros(len(Tsat_avg))
    rhol_avg = np.zeros(len(Tsat_avg))
    rhov_avg = np.zeros(len(Tsat_avg))
    rhor_avg = np.zeros(len(Tsat_avg))
    Zv_avg = np.zeros(len(Tsat_avg))
    Psat_95 = np.zeros(len(Tsat_avg))
    logPsat_95 = np.zeros(len(Tsat_avg))
    rhol_95 = np.zeros(len(Tsat_avg))
    rhov_95 = np.zeros(len(Tsat_avg))
    rhor_95 = np.zeros(len(Tsat_avg))
    Zv_95 = np.zeros(len(Tsat_avg))

    f0 = open('C:/Users/rmesserl/Documents/Helium_ab_initio/Critical_constants/'+system+'_GEMC.txt','w')
    f0.write('T (K)'+'\t'+'rhol (kg/m3)'+'\t'+'rhov (kg/m3)'+'\t'+'Psat (MPa)'+'\n')
    
    f1 = open('C:/Users/rmesserl/Documents/Helium_ab_initio/Critical_constants/'+system+'_GEMC_molar.txt','w')
#    f1.write('T (K)'+'\t'+'rhol (kmol/m3)'+'\t'+'rhov (kmol/m3)'+'\t'+'Psat (MPa)'+'\n')
    f1.write('T (K)'+'\t'+'rhol (kmol/m3)'+'\t'+'rhov (kmol/m3)'+'\t'+'Psat (MPa)'+'\t'+'Zvap'+'\n')
    
    for iTsat, Tsat_i in enumerate(Tsat_avg):
        Psat_avg[iTsat] = np.mean(Psat[Tsat == Tsat_i])
        rhol_avg[iTsat] = np.mean(rhol[Tsat == Tsat_i])
        rhov_avg[iTsat] = np.mean(rhov[Tsat == Tsat_i])
        rhor_avg[iTsat] = np.mean(rhor[Tsat == Tsat_i])
        Zv_avg[iTsat] = np.mean(Zv[Tsat == Tsat_i])
        
        Psat_95[iTsat] = 3.18 * np.std(Psat[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        logPsat_95[iTsat] = 3.18 * np.std(logPsat[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        rhol_95[iTsat] = 3.18 * np.std(rhol[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        rhov_95[iTsat] = 3.18 * np.std(rhov[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        rhor_95[iTsat] = 3.18 * np.std(rhor[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        Zv_95[iTsat] = 3.18 * np.std(Zv[Tsat == Tsat_i]) / len(Tsat[Tsat == Tsat_i])
        
        f0.write(str(Tsat_i)+'\t'+str(np.round(rhol_avg[iTsat],2))+'\t'+str(np.round(rhol_95[iTsat],2))+'\t'+str(np.round(rhov_avg[iTsat],2))+'\t'+str(np.round(rhov_95[iTsat],2))+'\t'+str(np.round(Psat_avg[iTsat]/10.,4))+'\t'+str(np.round(Psat_95[iTsat]/10.,4))+'\n')
#        f1.write(str(Tsat_i)+'\t'+str(np.round(rhol_avg[iTsat]/Mw_He,3))+'\t'+str(np.round(rhol_95[iTsat]/Mw_He,3))+'\t'+str(np.round(rhov_avg[iTsat]/Mw_He,3))+'\t'+str(np.round(rhov_95[iTsat]/Mw_He,3))+'\t'+str(np.round(Psat_avg[iTsat]/10.,4))+'\t'+str(np.round(Psat_95[iTsat]/10.,4))+'\n')
        f1.write(str(Tsat_i)+' & '+str(np.round(rhol_avg[iTsat]/Mw_He,3))+' $\pm$ '+str(np.round(rhol_95[iTsat]/Mw_He,3))+' & '+str(np.round(rhov_avg[iTsat]/Mw_He,4))+' $\pm$ '+str(np.round(rhov_95[iTsat]/Mw_He,4))+' & '+str(np.round(Psat_avg[iTsat]/10.,5))+' $\pm$ '+str(np.round(Psat_95[iTsat]/10.,5))+' & '+str(np.round(Zv_avg[iTsat],4))+' $\pm$ '+str(np.round(Zv_95[iTsat],4))+'\n')
        
    f0.close()    
    f1.close()
        
    logPsat_avg = np.log10(Psat_avg)
        
    He_fit = VLE(Tsat,rhol,rhov,Psat,Mw_He)
    
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
    Zc = He_fit.Zc
    
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
    
    urhoc, uTc, uPc, uZc = He_fit.bootstrapCriticals()
    ulogPc = (np.log10(Pc+uPc) - np.log10(Pc-uPc))/2.
    uinvTc = (10./(Tc-uTc)-10./(Tc+uTc))/2.
             
    Tplot_low = 10
    invTplot_low = 10./Tplot_low
    
    print('Critical temperature = '+str(np.round(Tc,2))+r'$ \pm$ '+str(np.round(uTc,2))+' K, Critical Pressure = '+str(np.round(Pc/10.,3))+r'$ \pm$ '+str(np.round(uPc/10.,3))+' MPa, Critical Density = '+str(np.round(rhoc/Mw_He,2))+' $\pm$ '+str(np.round(urhoc/Mw_He,2))+' (kmol/m$^3$), Critical Compressibility = '+str(np.round(Zc,4))+r'$ \pm$ '+str(np.round(uZc,4)))
    print('Critical temperature = '+str(np.round(Tc,2))+r'$ \pm$ '+str(np.round(uTc,2))+' K, Critical Pressure = '+str(np.round(Pc/10.,3))+r'$ \pm$ '+str(np.round(uPc/10.,3))+' MPa, Critical Density = '+str(np.round(rhoc,1))+' $\pm$ '+str(np.round(urhoc,2))+' (kg/m$^3$), Critical Compressibility = '+str(np.round(Zc,4))+r'$ \pm$ '+str(np.round(uZc,4)))
  
    ax0[0].errorbar(rhol_avg[Tsat_avg >= Tplot_low]/Mw_He,Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhol_95[Tsat_avg >= Tplot_low]/Mw_He,fmt=color+shape,mfc='None',markersize=10,label=system_label+', GEMC',capsize=6)
    ax0[0].errorbar(rhov_avg[Tsat_avg >= Tplot_low]/Mw_He,Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhov_95[Tsat_avg >= Tplot_low]/Mw_He,fmt=color+shape,mfc='None',markersize=10,capsize=6)
    ax0[0].errorbar(rhor_avg[Tsat_avg >= Tplot_low]/Mw_He,Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhor_95[Tsat_avg >= Tplot_low]/Mw_He,fmt=color+shape,mfc='None',markersize=10,capsize=6)
    ax0[0].plot(rhorplot[Tplot >= Tplot_low]/Mw_He,Tplot[Tplot >= Tplot_low],color+line,label=system_label+', Fit, GEMC')
    ax0[0].plot(rholRSplot[Tplot >= Tplot_low]/Mw_He,Tplot[Tplot >= Tplot_low],color+line)
    ax0[0].plot(rhovRSplot[Tplot>= Tplot_low]/Mw_He,Tplot[Tplot >= Tplot_low],color+line)
    ax0[0].errorbar(rhoc/Mw_He,Tc,xerr=urhoc/Mw_He,yerr=uTc,fmt=color+shape,markersize=10,mfc=color,label=system_label+', Critical, GEMC',capsize=6)
   
### Old, mass basis     
#    ax0[0].errorbar(rhol_avg[Tsat_avg >= Tplot_low],Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhol_95[Tsat_avg >= Tplot_low],fmt=color+shape,mfc='None',markersize=10,label=system_label+', GEMC',capsize=6)
#    ax0[0].errorbar(rhov_avg[Tsat_avg >= Tplot_low],Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhov_95[Tsat_avg >= Tplot_low],fmt=color+shape,mfc='None',markersize=10,capsize=6)
#    ax0[0].errorbar(rhor_avg[Tsat_avg >= Tplot_low],Tsat_avg[Tsat_avg >= Tplot_low],xerr=rhor_95[Tsat_avg >= Tplot_low],fmt=color+shape,mfc='None',markersize=10,capsize=6)
##    ax0[0].plot(rhol,Tsat,color+shape,mfc='None',label=system_label+', GEMC')
##    ax0[0].plot(rhov,Tsat,color+shape,mfc='None')
##    ax0[0].plot(rhor,Tsat,color+shape,mfc='None')
#    ax0[0].plot(rhorplot[Tplot >= Tplot_low],Tplot[Tplot >= Tplot_low],'k'+line,label=system_label+', Fit, GEMC')
#    ax0[0].plot(rholRSplot[Tplot >= Tplot_low],Tplot[Tplot >= Tplot_low],'k'+line)
#    ax0[0].plot(rhovRSplot[Tplot>= Tplot_low],Tplot[Tplot >= Tplot_low],'k'+line)
##    ax0[0].errorbar(rhoc,Tc,xerr=urhoc,yerr=uTc,fmt=color+'d',markersize=12,mfc='None',label=system_label+', Critical, GEMC')
#    ax0[0].errorbar(rhoc,Tc,xerr=urhoc,yerr=uTc,fmt=color+shape,markersize=10,mfc=color,label=system_label+', Critical, GEMC',capsize=6)
    
    ax0[1].errorbar(invTsat_avg[invTsat_avg <= invTplot_low],logPsat_avg[invTsat_avg <= invTplot_low]-1,yerr=logPsat_95[invTsat_avg <= invTplot_low],fmt=color+shape,mfc='None',markersize=10,capsize=6)#,label=system_label+', GEMC')
#    ax0[1].plot(invTsat,logPsat-1,color+shape,mfc='None',label=system_label+', GEMC')
    ax0[1].plot(invTplot[invTplot <= invTplot_low],logPsatplot[invTplot <= invTplot_low]-1,'k'+line)#,label=system_label+', Fit, GEMC')
#    ax0[1].errorbar(10./Tc,np.log10(Pc)-1,xerr=uinvTc,yerr=ulogPc,fmt=color+'d',markersize=14,mfc='None')#,label=system_label+', Critical, GEMC')
    ax0[1].errorbar(10./Tc,np.log10(Pc)-1,xerr=uinvTc,yerr=ulogPc,fmt=color+shape,markersize=10,mfc=color,capsize=6)#,label=system_label+', Critical, GEMC')
        
    ax0[1].plot([],[],color+shape,mfc='None',label=system_label+', GEMC',markersize=10)
#    ax0[1].plot([],[],color+'d',markersize=12,mfc='None',label=system_label+', Critical, GEMC')
    ax0[1].plot([],[],color+shape,markersize=12,mfc=color,label=system_label+', Critical, GEMC')
    
#    ax1.plot(Tsat,Zv,color+shape,mfc='None',markersize=10)
    ax1.errorbar(Tsat_avg,Zv_avg,yerr=Zv_95,fmt=color+shape,mfc='None',markersize=10,capsize=6)
    ax1.errorbar(Tc,Zc,xerr=uTc,yerr=uZc,fmt=color+shape,markersize=10,mfc=color,capsize=6)
    ax1.plot([],[],color+shape,mfc='None',label=system_label+', GEMC',markersize=10)
    ax1.plot([],[],color+shape,markersize=10,mfc=color,label=system_label+', Critical, GEMC')

ax0[0].errorbar(rhoc_Kofke/Mw_He,Tc_Kofke,xerr=urhoc_Kofke/Mw_He,yerr=uTc_Kofke,color='cyan',marker='d',mfc='cyan',markersize=10,label='Critical, VEOS',capsize=6)
#ax0[0].set_xlabel(r'$\rho$ (kmol/m$^3$)',fontsize=28)    
##ax0[0].errorbar(rhoc_Kofke,Tc_Kofke,xerr=urhoc_Kofke,yerr=uTc_Kofke,fmt='gv',markersize=10,mfc='g',label='Critical, VEOS',capsize=6)
##ax0[0].set_xlabel(r'$\rho$ (kg/m$^3$)',fontsize=28)
#ax0[0].set_ylabel('$T$ (K)',fontsize=28)

ax0[1].errorbar(10./Tc_Kofke,np.log10(Pc_Kofke)-1,xerr=uinvTc_Kofke,yerr=ulogPc_Kofke,color='cyan',marker='d',mfc='cyan',markersize=10,capsize=6)#,label='2-body, Critical, VEOS')
ax0[1].plot([],[],color='cyan',marker='d',linestyle='',mfc='cyan',markersize=10,label=r'(2+3)-body, Critical, VEOS$\infty$')
#ax0[1].set_xlabel(r'10/$T$ (1/K)',fontsize=28)
#ax0[1].set_ylabel(r'$\log_{10}$($P$/MPa)',fontsize=28)
ax0[1].legend(frameon=False)

#plt.tight_layout()

fig0.savefig('GEMC_results.pdf')
#plt.show()

ax1.errorbar(Tc_Kofke,Zc_Kofke,xerr=uTc_Kofke,yerr=uZc_Kofke,color='cyan',marker='d',markersize=10,mfc='cyan',capsize=6)
ax1.plot([],[],color='cyan',marker='d',linestyle='',markersize=10,mfc='cyan',label=r'(2+3)-body, Critical, VEOS$\infty$')

#ax1.set_ylabel(r'$Z_{\rm vap}^{\rm sat}$')
#ax1.set_xlabel(r'$T$ (K)')
ax1.legend(frameon=False)

#plt.tight_layout()

fig1.savefig('GEMC_Compressibility.pdf')
plt.show()