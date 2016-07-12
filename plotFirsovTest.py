"""
 plotFirsovTest.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the angular distribution from the Firsov slab test
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotFirsovTest.py filedir N_file N_in theta_0 energy_0 angle_bin
 ---------------------------------------------------------------------------------
 Parameters:
 - filedir = input path (string)
 - N_file = number of simulated files
 - N_in = number of simulated particles
 - theta_0 = incoming angle
 - energy_0 = incoming energy
 - angle_bin = bin angle
 --------------------------------------------------------------------------------
 Caveats:
 None
 ---------------------------------------------------------------------------------
 Modification history:
 - 2015/09/22: creation date
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
energy_0 = float(arg_list[5])
angle_bin = float(arg_list[6])


sphere_vol_id = 1002

# set-down
vecThetaOut = []
vecPhiOut = []
vecPsiOut = []
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

    # ---------------------------------


    where_sphere = np.where(vol_id == sphere_vol_id)
    evt_id_sphere = evt_id[where_sphere]
    ene_sphere = ene_ent[where_sphere]
    mdz_sphere = MDZ_ent[where_sphere]
    mdy_sphere = MDY_ent[where_sphere]
    mdx_sphere = MDX_ent[where_sphere]
    
    for jev_sphere in xrange(len(evt_id_sphere)):
		vecEnergyOut.append(ene_sphere[jev_sphere])
		vecThetaOut.append((180./np.pi)*np.arccos(-mdz_sphere[jev_sphere])-90.)
		vecPhiOut.append((180./np.pi)*np.arctan(mdy_sphere[jev_sphere]/mdx_sphere[jev_sphere]))
		vecPsiOut.append(((180./np.pi)*np.arccos(-mdz_sphere[jev_sphere])-90.)/theta_0)
	



# SPHERE
N_out = len(vecEnergyOut)
vecEnergyOut = np.array(vecEnergyOut)
vecPsiOut = np.array(vecPsiOut)

if (N_out > 1):
    angle_min = np.min(vecPsiOut)
    angle_max = np.max(vecPsiOut)
    n_bins_angle = int((angle_max - angle_min)/angle_bin)
    if (n_bins_angle < 1):
    	n_bins_angle = 100
else:
    print "NO PARTICLES TO PLOT!!!!"


def firsov_law(psi):
	return (3./(2.*np.pi))*((psi**(3./2.))/(1. + psi**3))
	
# histograms for the slab
N_array_out, bin_array_out = np.histogram(vecPsiOut, bins = n_bins_angle, range = (0, angle_max))
N_array_out_norm = np.zeros(len(N_array_out))
err_N_array_out = np.zeros(len(N_array_out))
err_angle = np.zeros(len(N_array_out))
angle_array = np.zeros(len(N_array_out))
left_angle_array = np.zeros(len(N_array_out))

N_model = np.zeros(len(N_array_out))


for jn in xrange(len(N_array_out)):
	err_angle[jn] = (bin_array_out[jn+1] - bin_array_out[jn])/2.
	N_array_out_norm[jn] = (float(N_array_out[jn])/(2.*err_angle[jn]))/N_in
	angle_array[jn] = bin_array_out[jn] + err_angle[jn]
	left_angle_array[jn] = bin_array_out[jn]
	err_N_array_out[jn] = (math.sqrt(float(N_array_out[jn]))/(2.*err_angle[jn]))/N_in
	N_model[jn] = firsov_law(angle_array[jn])


# Calculate degrees of freedom of fit
dof = len(angle_array)

# Plot the results
gs = gridspec.GridSpec(2, 1, height_ratios=[4,1]) 
gs.update(hspace=0.0)
ax = plt.subplot(gs[0])
ax1 = plt.subplot(gs[1])

# out angle distribution
ax.bar(left_angle_array, N_array_out_norm, width=2.*err_angle, edgecolor="lightblue", facecolor='lightblue', lw = 2, label='BoGEMMS')
ax.errorbar(angle_array, N_array_out_norm, xerr=err_angle, yerr=err_N_array_out, capsize=0, fmt='none', lw = 0.5, ecolor='black')
ax.set_xticklabels([])
ax.plot(angle_array, N_model, label='Firsov model')

ax.set_title('Soft proton slab test - E = '+str(energy_0)+' keV, $\Theta_{0}$ = '+str(theta_0))
ax1.set_xlabel("$\Psi$")
ax.set_ylabel("W($\Psi$)")
#ax.set_xlim(0, 1000)
#ax.set_ylim(np.min(N_array_slab)-np.min(N_array_slab)/10., np.max(N_array_slab)+np.max(N_array_slab)/10.)
ax.set_xlim(0, 30)
#ax.set_ylim(0, 1)
#ax.set_xscale('log')
#ax.set_yscale('log')

ax_range = ax.get_xlim()
ax1.set_xlim(ax_range)
ax1.errorbar(angle_array, ((N_array_out_norm-N_model)/N_array_out_norm)*100., xerr=err_angle, yerr=(err_N_array_out/N_array_out_norm)*100., fmt='.b', lw = 2, ecolor='black')
ax1.plot([0, ax_range[1]], [0,0], '-k')
ax1.set_ylabel("Residual [%]")
ax1.set_ylim(-100., 95.)

#ax.errorbar(angle_array, N_array_out, xerr=err_angle)

#ax.bar(left_ene_array_slab, N_array_slab, width=2.*err_ene_slab, yerr=err_N_array_slab, edgecolor="black", facecolor='white', lw = 2, ecolor='gray', label='Energy spectrum')
#plt.text(0.7, 0.94, 'Efficiency = ('+str(round(att_eff*100.,1))+' +/- '+str(round(err_att_eff*100.,3))+') %', transform=ax.transAxes, fontsize=12, zorder=100)


ax.legend(numpoints=1, loc=1)



plt.show()
hdulist.close()
