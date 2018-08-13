# -*- coding: utf-8 -*-
"""
Created on Wed Nov 15 14:31:30 2017

@author: ram9
"""

import numpy as np
import matplotlib.pyplot as plt
#from helium_2body_potential import V_total_hat

rcut = 14. #[A]

xyz = np.loadtxt('initial_configuration_debug_vapor5.txt',skiprows=1)
Lbox = 155.03813 #[A]
halfLbox = Lbox/2.

Nmol = len(xyz)

rmatrix = np.zeros([Nmol,Nmol])

for i in range(Nmol):
    
    for j in range(Nmol):
        
        dx = xyz[i,0] - xyz[j,0]
        dy = xyz[i,1] - xyz[j,1]
        dz = xyz[i,2] - xyz[j,2]
        
        if dx > halfLbox:
            
            dx = Lbox - dx
            
        if dy > halfLbox:
            
            dy = Lbox - dy
            
        if dz > halfLbox:
            
            dz = Lbox - dz
        
        dr = np.sqrt(dx**2 + dy**2 + dz**2)
        
        rmatrix[i,j] = dr
               
U_total = 0

for i in range(Nmol):
    
    for j in range(Nmol):
        
        if j > i:
            
            drij = rmatrix[i,j]
            
#            if drij < 1.1375:
#                
#                U_total += 452459.440452
#                
#            elif drij < 1.5:
#                
#                U_total += 452459.440452 + (5557.32389591 - 452459.440452)/(1.5-1.1375)*(drij-1.1375)
#                
#            elif drij < 1.8625:
#                
#                U_total += 5557.32389591 + (1057.18832413-5557.32389591)/(1.8625-1.5)*(drij-1.5)
#                
#            elif drij <= rcut:
#            
#                U_total += V_total_hat(rmatrix[i,j])
                
            if drij <= rcut:
                
                U_total += V_total_hat(rmatrix[i,j])
            
print(U_total)
print(U_total/Nmol)