"""
 testRemiz.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the angular distribution from the Remizovich slab test
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotRemizTest.py filedir N_file N_in theta_0 angle_bin
 ---------------------------------------------------------------------------------
 Parameters:
 - filedir = input path (string)
 - N_file = number of simulated files
 - N_in = number of simulated particles
 - theta_0 = incoming angle
 - angle_bin = bin angle
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2016/06/01: creation date
"""



import pyfits
import numpy as np
import math
import sys
import matplotlib.pyplot as plt
from matplotlib import gridspec

# Import the input parameters
arg_list = sys.argv
filedir = arg_list[1]
N_fits = int(arg_list[2])
N_in = int(arg_list[3])
theta_0 = float(arg_list[4])
angle_bin = float(arg_list[5])



def firsov_law(psi):
	return (3./(2.*np.pi))*((psi**(3./2.))/(1. + psi**3))
	
def remiz_law(psi, chi):
	Omega = np.sqrt((3.0*Psi)/((Psi**2.) - Psi + 1.0 + ((Chi/2.)**2.)))
	print "Omega", Omega
	return (1./(12.*(np.pi**2.)*(Psi**(1./2.))))*(((Omega**4.)/(1. + (Omega**2.))) + (Omega**3.)*np.arctan(Omega))
	
nbins = 100
#theta = np.linspace(10**(-15), 90., nbins)
#beta = np.linspace(-90., +90., nbins)
theta = np.linspace(10**(-15), np.pi/2., nbins)
beta = np.linspace(-np.pi/2., +np.pi/2., nbins)
theta_0 = theta_0*(np.pi/180.)
W = np.zeros(nbins*nbins)
for i in xrange(nbins):
	Psi = theta[i]/theta_0
	for j in xrange(nbins):
		ij = nbins*j + i
		Chi = beta[j]/theta_0
		print remiz_law(Psi, Chi)
		#W[ij] = remiz_law(Psi, Chi)
		#print W[ij]
