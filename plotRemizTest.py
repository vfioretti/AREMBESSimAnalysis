"""
 plotRemizTest.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the angular distribution from the Remizovich slab test
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotRemizTest.py filedir N_file N_in theta_0 energy_0 angle_bin
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
		vecPsiOut.append(((180./np.pi)*np.arccos(-mdz_sphere[jev_sphere])-90.)/theta_0)
		phi = (180./np.pi)*np.arctan(mdy_sphere[jev_sphere]/mdx_sphere[jev_sphere])
		vecPhiOut.append(phi)
		vecChiOut.append((phi)/theta_0)

	


# SPHERE
N_out = len(vecEnergyOut)
vecEnergyOut = np.array(vecEnergyOut)
vecPsiOut = np.array(vecPsiOut)
vecChiOut = np.array(vecChiOut)


if (N_out > 1):
    angle_min = np.min(vecPsiOut)
    angle_max = np.max(vecPsiOut)
    n_bins_psi = int((10.)/angle_bin)
    if (n_bins_psi < 1):
    	n_bins_psi = 100
    
    angle_min = np.min(vecChiOut)
    angle_max = np.max(vecChiOut)
    n_bins_chi = int((10.)/angle_bin)
    if (n_bins_chi < 1):
    	n_bins_chi = 100
else:
    print "NO PARTICLES TO PLOT!!!!"


def firsov_law(psi):
	return (3./(2.*np.pi))*((psi**(3./2.))/(1. + psi**3))
	
def remiz_law(Psi, Chi):
	Omega = np.sqrt((3.0*Psi)/((Psi**2.) - Psi + 1.0 + ((Chi/2.)**2.)))
	return (1./(12.*(np.pi**2.)*(Psi**(1./2.))))*(((Omega**4.)/(1. + (Omega**2.))) + (Omega**3.)*np.arctan(Omega))

print vecChiOut
print vecPsiOut
# histograms for the slab
N_histo2d = np.histogram2d(vecChiOut, vecPsiOut, bins = (n_bins_chi, n_bins_psi),  range=[[-5, 5.], [0, 10]])

N_array_out = N_histo2d[0]
extent_chi = N_histo2d[1]
extent_psi = N_histo2d[2]
bin_psi = extent_psi[1] - extent_psi[0]
bin_chi = extent_chi[1] - extent_chi[0]
pix_area = bin_chi*bin_psi


N_array_out_norm = N_array_out
for xel in xrange(len(extent_chi)-1):
	row = N_array_out[xel]
	temp_row_norm = []
	for yel in xrange(len(extent_psi)-1):
		temp_row_norm.append(((row[yel])/pix_area)/N_in)
	N_array_out_norm[xel] = temp_row_norm


N_array_model = np.zeros((len(extent_psi)-1, len(extent_chi)-1))
for yel in xrange(len(extent_psi)-1):
	for xel in xrange(len(extent_chi)-1):
		N_array_model[yel][xel] = remiz_law(extent_psi[yel]+bin_psi/2., extent_chi[xel]+bin_chi/2.)
	
print np.max(N_array_out_norm)
print np.max(N_array_model)
	
	
# Plotting

fig, axes = plt.subplots(1, 2, facecolor='w', figsize=(15,7), sharey=True)
fig.suptitle('Soft proton slab test - E = '+str(energy_0)+' keV, $\Theta_{0}$ = '+str(theta_0), fontsize=12)

ax_list = axes.flat

im1 = ax_list[0].imshow(N_array_out_norm.transpose()[::-1], extent = (extent_chi[0], extent_chi[len(extent_chi)-1], extent_psi[0], extent_psi[len(extent_psi)-1]), interpolation="nearest")
im2 = ax_list[1].imshow(N_array_model, origin='lower', extent = (extent_chi[0], extent_chi[len(extent_chi)-1], extent_psi[0], extent_psi[len(extent_psi)-1]), interpolation="nearest")


# cax is used to plot a colorbar for each subplot
div = make_axes_locatable(ax_list[0])
cax = div.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im1, cax=cax)
cbar.set_label('W($\Psi$, $\chi$)')

# cax is used to plot a colorbar for each subplot
div = make_axes_locatable(ax_list[1])
cax = div.append_axes("right", size="5%", pad=0.05)
cbar = plt.colorbar(im2, cax=cax)
cbar.set_label('W($\Psi$, $\chi$)')

ax_list[0].set_title('BoGEMMS')
ax_list[1].set_title('Remizovich Model')

ax_list[0].set_xlabel('$\chi$')
ax_list[0].set_ylabel('$\Psi$')
ax_list[1].set_xlabel('$\chi$')
ax_list[1].set_ylabel('$\Psi$')


#cax,kw = mpl.colorbar.make_axes([ax for ax in axes.flat])
#plt.colorbar(im, cax=cax, **kw)

#cbar1 = fig.colorbar(im, orientation='horizontal')

#cbar1 = fig.colorbar(im, orientation='horizontal')



plt.tight_layout()

# moves subplots down slightly to make room for the figure title
plt.subplots_adjust(top=0.9)
plt.show()


hdulist.close()
