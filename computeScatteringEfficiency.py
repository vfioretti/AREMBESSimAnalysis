"""
 plotScatteringEfficency2.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the soft proton scattering efficiency
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotScatteringEfficiency2.py filedir N_file N_in theta_0 energy_0 model
 ---------------------------------------------------------------------------------
 Parameters:
 - filedir = input path (string)
 - N_file = number of simulated files
 - N_in = number of simulated particles
 - theta_0 = incoming angle
 - energy_0 = incoming energy
 - angle_bin = bin angle
 - model = scattering model
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

model = arg_list[6]


solid_angle = 2.*np.pi*(1. - np.cos(np.arctan(0.6/933.)))
print "solid_angle = ", solid_angle

det_radius = (180./np.pi)*(np.arctan(0.6/933.))
angle_bin = det_radius*2.


target_vol_id = 1001
sphere_vol_id = 1002


N_in_eff = 0

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
    part_id = tbdata.field('PARTICLE_ID')
    
    # ---------------------------------


    where_target = np.where((vol_id == target_vol_id) & (part_id == 2212))
    evt_id_target = evt_id[where_target]
    ene_target = ene_ent[where_target]
    mdz_target = MDZ_ent[where_target]
    mdy_target = MDY_ent[where_target]
    mdx_target = MDX_ent[where_target]
    
    unique_evt_id_target = []
    jev_target = 0
    while (1):
		same_ev = np.where(evt_id_target == evt_id_target[jev_target])
		same_ev = same_ev[0]
		N_in_eff+= 1
		len_same_ev = len(same_ev)
		last_evt_id = same_ev[len_same_ev - 1]
		unique_evt_id_target.append(evt_id_target[jev_target])
		if (last_evt_id < (len(evt_id_target)-1)):
			jev_target = same_ev[len_same_ev - 1] + 1
		else:
			break
	
    
    where_sphere = np.where((vol_id == sphere_vol_id) & (ene_ent < energy_0) & (part_id == 2212))
    evt_id_sphere = evt_id[where_sphere]
    ene_sphere = ene_ent[where_sphere]
    mdz_sphere = MDZ_ent[where_sphere]
    mdy_sphere = MDY_ent[where_sphere]
    mdx_sphere = MDX_ent[where_sphere]
    
    for jev_sphere in xrange(len(evt_id_sphere)):
    	where_target = np.where(unique_evt_id_target == evt_id_sphere[jev_sphere])
    	if (where_target):
			vecEnergyOut.append(ene_sphere[jev_sphere])
			vecThetaOut.append((180./np.pi)*np.arccos(-mdz_sphere[jev_sphere])-90.)
			vecPsiOut.append(((180./np.pi)*np.arccos(-mdz_sphere[jev_sphere])-90.)/theta_0)
			phi = (180./np.pi)*np.arctan(mdy_sphere[jev_sphere]/mdx_sphere[jev_sphere])
			vecPhiOut.append(phi)
			vecChiOut.append((phi)/theta_0)

    		    		
# SPHERE
N_out = len(vecEnergyOut)
vecEnergyOut = np.array(vecEnergyOut)
vecPhiOut = np.array(vecPhiOut)
vecThetaOut = np.array(vecThetaOut)

angle_max = 5.

if (N_out > 1):
    angle_min_phi = np.min(vecPhiOut)
    angle_max_phi = angle_max
    n_bins_phi = int((angle_max - angle_min_phi)/angle_bin)
    if (n_bins_phi < 1):
    	n_bins_phi = 100
    
    angle_min_theta = np.min(vecThetaOut)
    angle_max_theta = angle_max
    n_bins_theta = int((angle_max - angle_min_theta)/angle_bin)
    if (n_bins_theta < 1):
    	n_bins_theta = 100
else:
    print "NO PARTICLES TO PLOT!!!!"


def firsov_law(psi):
	return (3./(2.*np.pi))*((psi**(3./2.))/(1. + psi**3))
	
def remiz_law(Psi, Chi):
	Omega = np.sqrt((3.0*Psi)/((Psi**2.) - Psi + 1.0 + ((Chi/2.)**2.)))
	return (1./(12.*(np.pi**2.)*(Psi**(1./2.))))*(((Omega**4.)/(1. + (Omega**2.))) + (Omega**3.)*np.arctan(Omega))


print "N_in: ", N_in
print "N_in_eff: ", N_in_eff

det_radius_rad = (det_radius*(np.pi/180.))
vecPhiOut_sel = vecPhiOut[np.where((vecPhiOut <= det_radius) & (vecPhiOut >= -det_radius))]
vecThetaOut_sel = vecThetaOut[np.where((vecPhiOut <= det_radius) & (vecPhiOut >= -det_radius))]


# Plot the results
gs = gridspec.GridSpec(1, 1, height_ratios=[4,1]) 
gs.update(hspace=0.0)
ax = plt.subplot(gs[0])



# Theta
N_array_out, bin_array_out = np.histogram(vecThetaOut_sel, bins = n_bins_theta, range = (0, angle_max))
N_array_out_norm = np.zeros(len(N_array_out))
err_N_array_out = np.zeros(len(N_array_out))
err_angle = np.zeros(len(N_array_out))
angle_array = np.zeros(len(N_array_out))
left_angle_array = np.zeros(len(N_array_out))
solid_angle_array = np.zeros(len(N_array_out))

for jn in xrange(len(N_array_out)):
	err_angle[jn] = (bin_array_out[jn+1] - bin_array_out[jn])/2.
	solid_angle_array[jn] = (2.*det_radius_rad)*(np.cos((90. - bin_array_out[jn+1])*(np.pi/180.)) - np.cos((90. - bin_array_out[jn])*(np.pi/180.)) )
	print solid_angle_array[jn]
	N_array_out_norm[jn] = (float(N_array_out[jn])/np.float(N_in_eff))/solid_angle_array[jn]
	angle_array[jn] = bin_array_out[jn] + err_angle[jn] + theta_0
	left_angle_array[jn] = bin_array_out[jn] + theta_0
	err_N_array_out[jn] = (np.sqrt(float(N_array_out[jn]))/np.float(N_in_eff))/solid_angle_array[jn]


# out angle distribution
ax.bar(left_angle_array, N_array_out_norm, width=2.*err_angle, edgecolor="lightblue", facecolor='lightblue', lw = 0.01, label='BoGEMMS', alpha=0.7)
ax.errorbar(angle_array, N_array_out_norm, xerr=err_angle, yerr=err_N_array_out, capsize=0, fmt='none', lw = 1, ecolor='black')

# experimental data
# 1MeV
if (energy_0 == 1000):
	if (theta_0 == 0.3):
		theta_exp = [0.63, 1.12, 1.61, 2.35, 3.08]
		eff_exp_max = [415.06 + 107.69, 584.64 + 151.61, 250.47 + 65.02, 98.88 + 25.71, 41.60 + 10.84]
		eff_exp_min = [415.06 - 107.69, 584.64 - 151.61, 250.47 - 65.02, 98.88 - 25.71, 41.60 - 10.84]
	if (theta_0 == 0.46):
		theta_exp = [0.63, 1.12, 1.61, 2.35, 3.08]
		eff_exp_max = [2.09 + 1.25, 413.98 + 107.47, 301.39 + 78.35, 116.20 + 30.21, 57.49 + 14.96]
		eff_exp_min = [2.09 - 1.25, 413.98 - 107.47, 301.39 - 78.35, 116.20 - 30.21, 57.49 - 14.96]
	if (theta_0 == 0.61):
		theta_exp = [1.12, 1.61, 2.35, 3.08]
		eff_exp_max = [160.94 + 41.89, 244.21 + 63.41, 136.86 + 35.61, 65.08 + 16.93]
		eff_exp_min = [160.94 - 41.89, 244.21 - 63.41, 136.86 - 35.61, 65.08 - 16.93]
	if (theta_0 == 0.89):
		theta_exp = [1.68, 2.42, 3.15, 4.13]
		eff_exp_max = [161.91 + 42.01, 114.94 + 29.82, 64.40 + 16.72, 30.44 + 7.93]
		eff_exp_min = [161.91 - 42.01, 114.94 - 29.82, 64.40 - 16.72, 30.44 - 7.93]
	if (theta_0 == 1.00):
		theta_exp = [2.35, 3.09, 4.07]
		eff_exp_max = [132.88 + 34.53, 80.00 + 20.80, 37.47 + 9.76]
		eff_exp_min = [132.88 - 34.53, 80.00 - 20.80, 37.47 - 9.76]
	if (theta_0 == 1.17):
		theta_exp = [2.35, 3.09, 4.07]
		eff_exp_max = [113.94 + 29.62, 84.25 + 21.91, 38.90 + 10.13]
		eff_exp_min = [113.94 - 29.62, 84.25 - 21.91, 38.90 - 10.13]



# Write output to file
name_fileout = "./diebold_test_results/"+str(int(energy_0))+"keV_"+str(theta_0)+"_"+str(N_in)+"/"+model+".dat"
print "Writing in "+name_fileout
f_out = open(name_fileout, 'wb')
for iel in xrange(len(angle_array)):
	# angle_x err_angle_x Eff err_Eff N_out N_in solid_angle 
	f_out.write(str(angle_array[iel])+" "+str(err_angle[iel])+" "+str(N_array_out_norm[iel])+" "+str(err_N_array_out[iel])+" "+str(N_array_out[iel])+" "+str(N_in_eff)+" "+str(solid_angle_array[iel])+"\n")


ax.fill_between(theta_exp, eff_exp_min, eff_exp_max)

ax.set_yscale("log")
ax.set_xlabel("Scattering angle [deg.]")
ax.set_ylabel("Scattering efficiency [sr$^{-1}$]")
ax.legend(numpoints=1, loc=1)
title = str(N_in)+" events, "+str(energy_0)+" keV, "+str(theta_0)+" deg., "+model+" model"
ax.set_title(title)

f_out.close()
plt.show()
hdulist.close()
