"""
 plotAngle.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the angular distribution from the Remizovich slab test
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotAngle.py filedir N_file N_in theta_0 energy_0 angle_bin model
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
import matplotlib as mpl
from mpl_toolkits.axes_grid1 import make_axes_locatable

# Import the input parameters
arg_list = sys.argv
filedir = arg_list[1]
N_fits = int(arg_list[2])
N_in = int(arg_list[3])
theta_0 = float(arg_list[4])
energy_0 = float(arg_list[5])
angle_bin = float(arg_list[6])
model = arg_list[7]


sphere_vol_id = 1002

# set-down
vecThetaOut = []
vecPhiOut = []
vecPsiOut = []
vecChiOut = []
vecEnergyOut = []


for jfits in xrange(N_fits):
    
    print '%%%%%%%%%%%%%% READING BoGEMMS FILE: ', jfits+1
    hdulist = pyfits.open(filedir+'/xyz.'+str(jfits)+'.fits.gz')
    
    tbdata = hdulist[1].data
    cols = hdulist[1].columns
    
    evt_id = tbdata.field('EVT_ID')
    vol_id = tbdata.field('VOLUME_ID')
    ene_ent = tbdata.field('E_KIN_ENT')
    MDZ_ent = tbdata.field('MDZ_ENT')
    MDY_ent = tbdata.field('MDY_ENT')
    MDX_ent = tbdata.field('MDX_ENT')
    X_ent = tbdata.field('X_ENT')
    part_id = tbdata.field('PARTICLE_ID')

    # ---------------------------------


    where_sphere = np.where((vol_id == sphere_vol_id) & (ene_ent < energy_0) & (part_id == 2212) & (X_ent < 0))
    evt_id_sphere = evt_id[where_sphere]
    ene_sphere = ene_ent[where_sphere]
    mdz_sphere = MDZ_ent[where_sphere]
    mdy_sphere = MDY_ent[where_sphere]
    mdx_sphere = MDX_ent[where_sphere]
    
    for jev_sphere in xrange(len(evt_id_sphere)):
		vecEnergyOut.append(ene_sphere[jev_sphere])
		vecThetaOut.append((180./np.pi)*np.arccos(-mdz_sphere[jev_sphere])-90.)
		vecPsiOut.append(((180./np.pi)*np.arccos(-mdz_sphere[jev_sphere])-90.)/theta_0)
		phi = (180./np.pi)*np.arctan(mdy_sphere[jev_sphere]/mdx_sphere[jev_sphere])
		vecPhiOut.append(phi)
		vecChiOut.append((phi)/theta_0)

	


# SPHERE
N_out = len(vecEnergyOut)
vecEnergyOut = np.array(vecEnergyOut)
vecThetaOut = np.array(vecThetaOut)
vecPhiOut = np.array(vecPhiOut)


if (N_out > 1):
    angle_min_phi = np.min(vecPhiOut)
    angle_max_phi = np.max(vecPhiOut)
    n_bins_phi = int((angle_max_phi - angle_min_phi)/angle_bin)
    if (n_bins_phi < 1):
    	n_bins_phi = 100
    
    angle_min_theta = np.min(vecThetaOut)
    angle_max_theta = np.max(vecThetaOut)
    n_bins_theta = int((angle_max_theta - angle_min_theta)/angle_bin)
    if (n_bins_theta < 1):
    	n_bins_theta = 100
else:
    print "NO PARTICLES TO PLOT!!!!"


def firsov_law(psi):
	return (3./(2.*np.pi))*((psi**(3./2.))/(1. + psi**3))
	
def remiz_law(Psi, Chi):
	Omega = np.sqrt((3.0*Psi)/((Psi**2.) - Psi + 1.0 + ((Chi/2.)**2.)))
	return (1./(12.*(np.pi**2.)*(Psi**(1./2.))))*(((Omega**4.)/(1. + (Omega**2.))) + (Omega**3.)*np.arctan(Omega))


# Plot the results
gs = gridspec.GridSpec(1, 1, height_ratios=[4,1]) 
gs.update(hspace=0.0)
ax = plt.subplot(gs[0])

# Phi
N_array_out, bin_array_out = np.histogram(vecPhiOut, bins = n_bins_phi, range = (0, angle_max_phi))
N_array_out_norm = np.zeros(len(N_array_out))
err_N_array_out = np.zeros(len(N_array_out))
err_angle = np.zeros(len(N_array_out))
angle_array = np.zeros(len(N_array_out))
left_angle_array = np.zeros(len(N_array_out))


for jn in xrange(len(N_array_out)):
	err_angle[jn] = (bin_array_out[jn+1] - bin_array_out[jn])/2.
	N_array_out_norm[jn] = (float(N_array_out[jn])/(2.*err_angle[jn]))/N_in
	angle_array[jn] = bin_array_out[jn] + err_angle[jn]
	left_angle_array[jn] = bin_array_out[jn]
	err_N_array_out[jn] = (math.sqrt(float(N_array_out[jn]))/(2.*err_angle[jn]))/N_in



# out angle distribution
ax.bar(left_angle_array, N_array_out_norm, width=2.*err_angle, edgecolor="lightblue", facecolor='lightblue', lw = 2, label='Phi', alpha=0.7)
ax.errorbar(angle_array, N_array_out_norm, xerr=err_angle, yerr=err_N_array_out, capsize=0, fmt='o', color='lightblue', lw = 1, ecolor='black')


# Theta
N_array_out, bin_array_out = np.histogram(vecThetaOut, bins = n_bins_theta, range = (0, angle_max_theta))
N_array_out_norm = np.zeros(len(N_array_out))
err_N_array_out = np.zeros(len(N_array_out))
err_angle = np.zeros(len(N_array_out))
angle_array = np.zeros(len(N_array_out))
left_angle_array = np.zeros(len(N_array_out))


for jn in xrange(len(N_array_out)):
	err_angle[jn] = (bin_array_out[jn+1] - bin_array_out[jn])/2.
	N_array_out_norm[jn] = (float(N_array_out[jn])/(2.*err_angle[jn]))/N_in
	angle_array[jn] = bin_array_out[jn] + err_angle[jn]
	left_angle_array[jn] = bin_array_out[jn]
	err_N_array_out[jn] = (math.sqrt(float(N_array_out[jn]))/(2.*err_angle[jn]))/N_in


# out angle distribution
ax.bar(left_angle_array, N_array_out_norm, width=2.*err_angle, edgecolor="yellow", facecolor='yellow', lw = 2, label='Theta', alpha=0.3)
ax.errorbar(angle_array, N_array_out_norm, xerr=err_angle, yerr=err_N_array_out, capsize=0, fmt='o', color='yellow', lw = 1, ecolor='black')

ax.set_xlabel("[deg]")
ax.set_ylabel("Norm. counts deg$^{-1}$")

ax.legend(numpoints=1, loc=1)
title = str(N_in)+" events, "+str(energy_0)+" keV, "+str(theta_0)+" deg., "+model+" model"
ax.set_title(title)

plt.show()
hdulist.close()
