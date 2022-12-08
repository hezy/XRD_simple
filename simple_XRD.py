#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Dec  7 14:56:13 2022

@author: hezy
"""

import numpy as np
import matplotlib.pyplot as plt
from scipy.special import wofz
from scipy.optimize import curve_fit
from itertools import combinations_with_replacement


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


def peaks_width(theta, U, V, W):
    theta_rad = theta*np.pi/180
    return np.sqrt( U*np.tan(theta_rad/2)**2 + V*np.tan(theta_rad/2) + W)


def intensity(theta_space, peaks_positions, peaks_width):
    y = np.zeros(500)
    for n in range(9):
        #print(n, peaks_position[n], peaks_width[n])
        y = y + peak(theta_space, peaks_positions[n], 1, peaks_width[n], 0.5)
    return y


def find_d(indices_list, a):
    miller = np.array(indices_list)
    return a/np.sqrt(miller.T[0]**2 + miller.T[1]**2 + miller.T[2]**2)


def bragg_angels(wavelength, d_spacings):
    return 2 * 180/np.pi * np.arcsin(wavelength/(2*d_spacings))  # *2 for 2θ  
    

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
    

N = 500
theta_space = np.linspace (0, 50, N)
#peaks_position = np.array([10.0,20.0,30.0,40.0,50.0,60.0,70.0,80.0,90.0])

wavelength = 0.15418  # CuKα radiation in nm
U, V, W = 0.2, 0.1, 0.05

#angular_intensity = intensity(theta_space, peaks_position, peaks_width(peaks_position, U, V, W))
#print(angular_intensity)

#make_graph(theta_space,angular_intensity)


### In simple cubic lattince, all Miller indices are allowed
sample_list = [0, 1, 2, 3]
SC_indices = list(combinations_with_replacement(sample_list, 3))
SC_indices.remove((0,0,0))
# print(SC_indices)
d_SC = find_d(SC_indices, 1)
# print(d_SC)

a_SC = 0.3352 # a for SC Polonium (α-Po), from https://en.wikipedia.org/wiki/Polonium
bragg_angels_SC = bragg_angels(wavelength, d_SC)

angular_intensity_SC = intensity(theta_space, bragg_angels_SC, peaks_width(bragg_angels_SC, U, V, W))

make_graph(theta_space,angular_intensity_SC)


### In body centerd cubic lattice, only indices with h+k+l=even are allowed
BCC_indices = SC_indices[:]
for item in BCC_indices:
        if (item[0] + item[1] + item[2]) % 2 != 0:
            BCC_indices.remove(item)
# print(BCC_indices)  
d_BCC = find_d(BCC_indices, 1)
# print(d_BCC)


### In face centered cubic lattice, h,k,l must all be either odd or even
FCC_indices = SC_indices[:]
for item in FCC_indices:
        all = "mixed"
        if (item[0]%2 != 0) and (item[1]%2 != 0) and (item[2]%2 != 0):
            all = "all pair"
        if (item[0]%2 == 0) and (item[1]%2 == 0) and (item[2]%2 == 0):
            all = "all even"
        if all == "mixed":
            FCC_indices.remove(item)
# print(FCC_indices)      
d_FCC = find_d(FCC_indices, 1)
# print(d_FCC)

