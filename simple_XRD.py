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


def lorentz(x, FWHM):
    """
    Lorentz function with max=1
    """  
    gamma = FWHM
    return 1 / (1 + np.square(x/gamma)) 


def gauss(x, FWHM):
    """
    Gauss function with max=1
    """
    sigma = FWHM/np.sqrt(2*np.log(2))
    return np.exp(- x**2 / (2* sigma**2))


def voigt(x, FWHM_L, FWHM_G):
    """
    Voigt function with max=1
    """
    gamma = FWHM_L
    sigma = FWHM_G/np.sqrt(2*np.log(2))
    z = (x + 1j*gamma)/np.sqrt(2)/sigma
    return np.sqrt(2*np.pi) * np.real(wofz(z))/np.sqrt(2*np.pi)/sigma
    # normolized Voigt (integral = 1): c * np.real(wofz((x + 1j*gamma)/(sigma * np.sqrt(2)))) / (sigma * np.sqrt(2*np.pi))
    # for Lorentz sigma=0, gamma=1, c=1
    # for Gauss sigma=1, gamma=0, c=1
  
    
def pseudo_voigt(x, FWHM, n):
    """
    pseudo-Voigt with max=1
    """
    return n * gauss(x, FWHM) + (1-n) * lorentz(x, FWHM)


def peak(x, x0, A, FWHM, n):
    """
    A single pseudo Voigt peak at x0, with hight A, width FWHM, and mixing factor n
    """
    return A * pseudo_voigt(x-x0, FWHM, n)


def peaks_width(two_theta, U, V, W):
    """
    Calaculate the width of a peak at a given theta
    """
    two_theta_rad = two_theta * np.pi/180
    return np.sqrt( U*np.tan(two_theta_rad/2)**2 + V*np.tan(two_theta_rad/2) + W)


def intensity(theta_space, peaks_positions, peaks_width, mixing_factor):
    """
    Adding up the contributions of all the peaks
    """
    y = np.zeros(len(theta_space))
    for n in range(len(peaks_positions)):
        #print(n, peaks_position[n], peaks_width[n])
        y = y + peak(theta_space, peaks_positions[n], 1, peaks_width[n], mixing_factor)
    return y


def find_d(indices_list, a):
    """
    Calculate plane distances for the various Miller indices
    """
    miller = np.array(indices_list)
    return a/np.sqrt(miller.T[0]**2 + miller.T[1]**2 + miller.T[2]**2)


def bragg_angels(wavelength, d_spacings):
    """
    Calculate the Bragg angles from the plane d spacings
    """
    return 2 * 180/np.pi * np.arcsin(wavelength/(2*d_spacings))  # *2 for 2θ  
    

def make_graph (x, y):
    """
    Plot a graph of intensity y vs angle x
    """
    fig1, ax = plt.subplots(figsize=(14, 8))
    ax.grid(visible=True, which='both', axis='both')
    ax.minorticks_on()
    ax.set_title("XRD", fontsize=16)
    ax.set_xlabel(r"$2 \theta$", fontsize=14)
    #ax.set_xlim()
    ax.set_ylabel(r"Intensity", fontsize=14)
    #ax.set_ylim()
    #ax.plot(x,y, '.r', label='experiment')
    ax.plot(x,y, '-b', label='theory')
    ax.legend()
    


N = 500
theta_space = np.linspace (0, 90, N)

wavelength = 0.15418  # CuKα radiation in nm
U, V, W = 0.1, 0.05, 0.02
sample_list = [0, 1, 2, 3, 4, 5, 6]  # for Miller indices


"""
SC - simple cubic lattice
all Miller indices are allowed
"""
# finding the Miller indices and planes distances for SC
SC_indices = list(combinations_with_replacement(sample_list, 3))
SC_indices.remove((0,0,0))
# print(SC_indices)
d_SC = find_d(SC_indices, 1)
# print(d_SC)

""" 
Particular case of SC - Polonium (α-Po)
data from https://en.wikipedia.org/wiki/Polonium
"""
a_SC = 0.3352 # (nm)
bragg_angels_SC = bragg_angels(wavelength, d_SC)
angular_intensity_SC = intensity(theta_space, bragg_angels_SC, peaks_width(bragg_angels_SC, U, V, W), 0.5)
make_graph(theta_space,angular_intensity_SC)


"""
BCC - body centerd cubic lattice
only indices with h+k+l=even are allowed
"""
BCC_indices = SC_indices[:]
for item in BCC_indices:
        if (item[0] + item[1] + item[2]) % 2 != 0:
            BCC_indices.remove(item)
# print(BCC_indices)  
d_BCC = find_d(BCC_indices, 1)
# print(d_BCC)

"""
Particular case of BCC - Tantalum (α-Ta)
Data from https://en.wikipedia.org/wiki/Tantalum
"""
a_BCC = 0.33058 #(nm)
bragg_angels_BCC = bragg_angels(wavelength, d_BCC)
angular_intensity_BCC = intensity(theta_space, bragg_angels_BCC, peaks_width(bragg_angels_BCC, U, V, W), 0.5)
make_graph(theta_space,angular_intensity_BCC)


"""
FCC - face centered cubic lattice
h,k,l must all be either odd or even
"""
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

"""
Particular case of FCC - Platinum
data from https://en.wikipedia.org/wiki/Lattice_constant
"""
a_FCC = 0.392 # (nm)
bragg_angels_FCC = bragg_angels(wavelength, d_FCC)
angular_intensity_FCC = intensity(theta_space, bragg_angels_FCC, peaks_width(bragg_angels_FCC, U, V, W), 0.5)
make_graph(theta_space,angular_intensity_FCC)

