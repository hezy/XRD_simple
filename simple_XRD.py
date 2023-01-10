#!/usr/bin/env python3
# -*- coding: utf-8 -*-

"""
Created on Wed Dec  7 14:56:13 2022
@author: Yehezkel Amiel

How is it now? (2023-01-10)
    seems to work fine for cubic symmetries (sc, bcc, fcc) with one element

Things to do (short term):
    [] function "intensity" should not use built in values
    [] function "make_graph" should give nicer titles
    [] add random noise
    [] add bacground
    [] investigate U,V,W values
    [] export to csv
    [] read setup and sample data from external files

For the future (long term):
    [] find Miller indecies based on lattice vectors rather than selection rules
    [] other symmetries
    [] more than one element    
"""

''' Importing lipreries '''
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz
#from scipy.optimize import curve_fit
#from itertools import combinations_with_replacement
from itertools import product


''' Functions '''

def lorentz(x, wL):
    # Lorentz with max=1 and w=FWHM: 
    gamma = wL
    return 1 / (1 + np.square(x/gamma)) 


def gauss(x, wG):
    # Gauss with max=1 and w=FWH
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


def peaks_width(theta, U, V, W):
    theta_rad = theta*np.pi/180
    return np.sqrt( U*np.tan(theta_rad/2)**2 + V*np.tan(theta_rad/2) + W)


def intensity(theta_space, peaks_positions, peaks_width):
    y = np.zeros(2000)
    for n in range(np.size(peaks_positions)):
        #print(n, peaks_positions[n], peaks_width[n])
        y = y + peak(theta_space, peaks_positions[n], 1, peaks_width[n], 0.5)
    return y


def find_d(indices_list, a):
    miller = np.array(indices_list)
    return a/np.sqrt(miller.T[0]**2 + miller.T[1]**2 + miller.T[2]**2)


def bragg_angels(wavelength, d_spacings):
    sintheta = wavelength/(2*d_spacings)
    sintheta = sintheta[abs(sintheta)<=1]  # removing values outside (-1,1)
    return 2 * 180/np.pi * np.arcsin(sintheta)  # *2 for 2θ  
    

def make_graph (x, y):
    fig1, ax = plt.subplots(figsize=(14, 8))
    ax.grid(visible=True, which='both', axis='both')
    ax.minorticks_on()
    ax.set_title("XRD", fontsize=16)
    ax.set_xlabel(r"$2 \theta$", fontsize=14)
    #ax.set_xlim()
    ax.set_ylabel(r"Intensity", fontsize=14)
    #ax.set_ylim()
    ax.plot(x,y, '.r', label='experiment')
    ax.plot(x,y, '-b', label='theory')
    ax.legend()
    
'''
===== 
Setup 
=====
'''

N = 2000
theta_space = np.linspace (0, 180, N)

wavelength = 0.15418  # CuKα radiation in nm
#wavelength = 0.18125  # 
U, V, W = 0.2, 0.2, 0.2


'''
============
Simple Cubic
============
'''

''' In simple cubic lattince, all Miller indices are allowed '''
sample_list = [-5,-4,-3,-2,-1,0,1,2,3,4,5]
indices_SC = list(product(sample_list, repeat = 3))
indices_SC.remove((0,0,0))


'''
Lattice parameter for SC Polonium (α-Po)
from https://en.wikipedia.org/wiki/Polonium 
'''
a_SC = 0.3352 

d_SC = find_d(indices_SC, a_SC)

bragg_angels_SC = bragg_angels(wavelength, d_SC)
angular_intensity_SC = intensity(theta_space,
                                 bragg_angels_SC,
                                 peaks_width(bragg_angels_SC, U, V, W))
make_graph(theta_space, angular_intensity_SC)


'''
===================
Body Centered Cubic
===================
'''

''' In body centerd cubic lattice, only indices with h+k+l=even are allowed '''
indices_BCC = []
for item in indices_SC:
        if (item[0] + item[1] + item[2]) % 2 == 0:
            indices_BCC.append(item)
# print(indices_BCC)  
 
'''
Lattice parameter for BCC Tantalum (α-Ta)
from https://en.wikipedia.org/wiki/Tantalum
'''
a_BCC = 0.33058

'''
Lattice parameter for BCC Tungsten (W)
from https://en.wikipedia.org/wiki/Lattice_constant
'''
#a_BCC = 0.3155

d_BCC = find_d(indices_BCC, a_BCC)
# print(d_BCC)

bragg_angels_BCC = bragg_angels(wavelength, d_BCC)
angular_intensity_BCC = intensity(theta_space,
                                  bragg_angels_BCC,
                                  peaks_width(bragg_angels_BCC, U, V, W))
make_graph(theta_space,angular_intensity_BCC)


'''
===================
Face Centered Cubic
===================
'''

''' In face centered cubic lattice, h,k,l must all be either odd or even '''
indices_FCC = []
#print('before:' ,indices_SC)
for item in indices_SC:
        if [(-1)**item[0], (-1)**item[1] ,(-1)**item[2]] == [1,1,1]:
            indices_FCC.append(item)
        if [(-1)**item[0], (-1)**item[1] ,(-1)**item[2]] == [-1,-1,-1]:
            indices_FCC.append(item)  

'''
Lattice parameter for FCC Platinum
from https://periodictable.com/Elements/078/data.html
'''
a_FCC = 0.39242 

'''
Lattice parameter a for FCC Pb
from https://en.wikipedia.org/wiki/Lattice_constant
'''
#a_FCC = 0.4920 

d_FCC = find_d(indices_FCC, a_FCC)
# print(d_FCC)

bragg_angels_FCC = bragg_angels(wavelength, d_FCC)
angular_intensity_FCC = intensity(theta_space,
                                  bragg_angels_FCC,
                                  peaks_width(bragg_angels_FCC, U, V, W))
make_graph(theta_space,angular_intensity_FCC)
