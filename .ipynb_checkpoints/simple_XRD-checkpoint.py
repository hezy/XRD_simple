#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:56:13 2022

@author: hezy
"""

import numpy as np
import matplotlib.pyplot as plt


def lorentz(x, wL):
    # Lorentz with max=1 and w=FWHM: 
    gamma = wL
    return 1 / (1 + np.square(x/gamma)) 


def gauss(x, wG):
    # Gauss with max=1 and w=FWHM
    sigma = wG/np.sqrt(2*np.log(2))
    return np.exp(- x**2 / (2* sigma**2))


def voigt(x, wL, wG):
    gamma = wL
    sigma = wG/np.sqrt(2*np.log(2))
    z = (x + 1j*gamma)/np.sqrt(2)/sigma
    return np.sqrt(2*np.pi) * np.real(wofz(z))/np.sqrt(2*np.pi)/sigma
    # normolized Voigt (integral = 1): c * np.real(wofz((x + 1j*gamma)/(sigma * np.sqrt(2)))) / (sigma * np.sqrt(2*np.pi))
    # for Lorentz sigma=0, gamma=1, c=1
    # for Gauss sigma=1, gamma=0, c=1
  
    
def pseudo_voigt(x, w, n):
    # pseudo-voigt with max=1 and w=FWHM:
    return n * gauss(x, w) + (1-n) * lorentz(x,w)


def peak(x, x0, A, w, n):
    return A * pseudo_voigt(x-x0, w, n)


xi = np.array([10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0])
xi_rad = xi*np.pi/180
U, V, W = 0.2, 0.1, 0.05
wi = np.sqrt( U*np.tan(xi_rad/2)**2 + V*np.tan(xi_rad/2) + W)

x = np.linspace (0, 120 , 500)
y = np.zeros(500)


for n in range(9):
    print(n, xi[n], wi[n])
    y = y + peak(x, xi[n], 1, wi[n], 0.5)

fig1, ax = plt.subplots(figsize=(14, 8))
ax.grid(visible=True, which='both', axis='both')
ax.minorticks_on()
ax.set_title("Gaussian and Lorentzian", fontsize=16)
ax.set_xlabel("x", fontsize=14)
#ax.set_xlim()
ax.set_ylabel("y", fontsize=14)
#ax.set_ylim()
ax.plot(x,y, '.r', label='experiment')
ax.plot(x,y, '-b', label='theory')
#ax.plot(x,yL, '-r', label='Lorentz')
#ax.plot(x,yG, '-b', label='Gauss')
#ax.plot(x,yPV, '-m', label='Pseudo Voigt')
#ax.plot(x,yV, '-g', label='Voigt')
ax.legend()


# In simple cubic lattince, all Miller indices are allowed
from itertools import combinations_with_replacement
sample_list = [0, 1, 2, 3]
SC_indices = list(combinations_with_replacement(sample_list, 3))
SC_indices.remove((0,0,0))
print(SC_indices)

def find_d(indices_list, a):
    miller = np.array(indices_list)
    return a/np.sqrt(miller.T[0]**2 + miller.T[1]**2 + miller.T[2]**2)

d_SC = find_d(SC_indices, 1)
print(d_SC)


# In body centerd cubic lattice, only indices with h+k+l=even are allowed
BCC_indices = SC_indices[:]
for item in BCC_indices:
        if (item[0] + item[1] + item[2]) % 2 != 0:
            BCC_indices.remove(item)
print(BCC_indices)  

d_BCC = find_d(BCC_indices, 1)
print(d_BCC)


#In face centered cubic lattice, h,k,l must all be either odd or even
FCC_indices = SC_indices[:]
for item in FCC_indices:
        all = "mixed"
        if (item[0]%2 != 0) and (item[1]%2 != 0) and (item[2]%2 != 0):
            all = "all pair"
        if (item[0]%2 == 0) and (item[1]%2 == 0) and (item[2]%2 == 0):
            all = "all even"
        if all == "mixed":
            FCC_indices.remove(item)
print(FCC_indices)    

d_FCC = find_d(FCC_indices, 1)
print(d_FCC)
