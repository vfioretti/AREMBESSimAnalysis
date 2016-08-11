"""
 plotTotalEnergyLoss.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the soft proton scattering efficiency
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotTotalEnergyLoss.py N_in energy_0 rebin model
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
N_in = int(arg_list[1])
energy_0 = float(arg_list[2])
rebin = int(arg_list[3])
model = arg_list[4]

# experimental data
# 1MeV
N_exp = 6
color_exp = ['Gray', 'red', 'green', 'blue', 'purple', 'OrangeRed']
if energy_0 == 1000:
    theta_0_array = [0.3, 0.46, 0.61, 0.89, 1.0, 1.17]
if energy_0 == 500:
    theta_0_array = [0.33, 0.48, 0.64, 0.85, 1.02, 1.19]
if energy_0 == 250:
    theta_0_array = [0.36, 0.51, 0.67, 0.89, 1.06, 1.23]


det_radius = (180./np.pi)*(np.arctan(0.6/933.))*1.
angle_bin = det_radius*2.
angle_max = 5.

handles2 = []


for jexp in xrange(N_exp):
    
    # Plot the results
    fig = plt.figure(1, figsize=(15, 10))
    ax2 = fig.add_subplot(111)

    # Reading simulation output
    name_file_in = './diebold_test_results/'+str(int(energy_0))+'keV_'+str(theta_0_array[jexp])+'_'+str(N_in)+'/'+model+'_energy.dat'
    f_in = open(name_file_in, 'r')
	
    energy_array_sim = []
    theta_out = []
    phi_out = []
    
    for line in f_in:
        line = line.strip()
        columns = line.split()
        # angle_x err_angle_x Eff err_Eff N_out N_in solid_angle
        columns[0] = float(columns[0])
        columns[1] = float(columns[1])
        columns[2] = float(columns[2])
        
        theta_out.append(columns[0])
        phi_out.append(columns[1])
        energy_array_sim.append(columns[2])
    
    if (len(energy_array_sim) > 1):  
    
       angle_min_theta = np.min(theta_out)
       angle_max_theta = angle_max
       n_bins_theta = int((angle_max - angle_min_theta)/angle_bin)
       if (n_bins_theta < 1):
    	   n_bins_theta = 100

    else:
        print "NO PARTICLES TO PLOT!!!!"

    N_theta_out, bin_theta_out = np.histogram(theta_out, bins = n_bins_theta, range = (0, angle_max))
    
    mean_energy_loss = []
    err_mean_energy_loss = []
    theta_xrange = []
    err_theta_xrange = []

    mean_energy_loss_rebin = []
    err_mean_energy_loss_rebin = []
    theta_xrange_rebin = []
    err_theta_xrange_rebin = []

    energy_array_sim = np.array(energy_array_sim)
    for jbin in xrange(n_bins_theta):
         theta_xrange.append(bin_theta_out[jbin] + (bin_theta_out[jbin+1] - bin_theta_out[jbin])/2.  + theta_0_array[jexp])
         err_theta_xrange.append((bin_theta_out[jbin+1] - bin_theta_out[jbin])/2.)
         theta_out = np.array(theta_out)
         energy_sel = energy_array_sim[np.where((theta_out > bin_theta_out[jbin]) & (theta_out < bin_theta_out[jbin+1]))]
         mean_energy_loss.append(energy_0 - np.mean(energy_sel))
         err_mean_energy_loss.append(np.std(energy_sel))

    jel=0
    while(1):
        if (jel + rebin) < len(theta_xrange):
            eloss = 0
            err_eloss_array = []
            err_theta_tot = 0
            for i in xrange(rebin):
                eloss+=mean_energy_loss[jel + i]
                err_eloss_array.append((err_mean_energy_loss[jel + i])**2.)
                err_theta_tot+=err_theta_xrange[jel + i]
            err_mean_energy_loss_rebin.append(np.sqrt(np.sum(err_eloss_array)))
            mean_energy_loss_rebin.append(eloss/np.float(rebin))
            theta_xrange_rebin.append(theta_xrange[jel] + (theta_xrange[jel+ rebin-1] - theta_xrange[jel])/2.)
            err_theta_xrange_rebin.append(err_theta_tot)
        else:
            break
        jel = jel + rebin


    if (energy_0 == 1000):
	if (theta_0_array[jexp] == 0.3):
		theta_exp = [0.63, 1.12, 1.61, 2.35, 3.08]
		ene_loss_max = [32. + 2, 39. + 2., 40. + 2., 47. + 2., 53. + 2.]
		ene_loss_min = [32. - 2, 39. - 2., 40. - 2., 47. - 2., 53. - 2.]
	if (theta_0_array[jexp] == 0.46):
		theta_exp = [1.12, 1.61, 2.35, 3.08]
		ene_loss_max = [39. + 2., 40. + 2., 47. + 2., 53. + 2.]
		ene_loss_min = [39. - 2., 40. - 2., 47. - 2., 53. - 2.]
	if (theta_0_array[jexp] == 0.61):
		theta_exp = [1.12, 1.61, 2.35, 3.08]
		ene_loss_max = [38. + 2., 39. + 2., 45. + 2., 49. + 2.]
		ene_loss_min = [38. - 2., 39. - 2., 45. - 2., 49. - 2.]
	if (theta_0_array[jexp] == 0.89):
		theta_exp = [1.68, 2.42, 3.15, 4.13]
		ene_loss_max = [38. + 2., 39. + 2., 44. + 2., 49. + 2.]
		ene_loss_min = [38. - 2., 39. - 2., 44. - 2., 49. - 2.]
	if (theta_0_array[jexp] == 1.00):
		theta_exp = [2.35, 3.09, 4.07]
		ene_loss_max = [39. + 2., 44. + 2., 52. + 2.]
		ene_loss_min = [39. - 2., 44. - 2., 52. - 2.]
	if (theta_0_array[jexp] == 1.17):
		theta_exp = [2.35, 3.09, 4.07]
		ene_loss_max = [39. + 2., 46. + 2., 54. + 2.]
		ene_loss_min = [39. - 2., 46. - 2., 54. - 2.]
    if (energy_0 == 250):
	if (theta_0_array[jexp] == 0.36):
		theta_exp = [0.69, 1.18, 1.67, 2.41, 3.14]
                ene_loss_max = [13. + 10, 22. + 10., 32. + 10., 35. + 10., 40. + 10.]
                ene_loss_min = [13. - 10, 22. - 10., 32. - 10., 35. - 10., 40. - 10.]
	if (theta_0_array[jexp] == 0.51):
		theta_exp = [1.18, 1.67, 2.41, 3.14]
                ene_loss_max = [22. + 10, 30. + 10., 38. + 10., 36. + 10.]
                ene_loss_min = [22. - 10, 30. - 10., 38. - 10., 36. - 10.]
	if (theta_0_array[jexp] == 0.67):
		theta_exp = [1.18, 1.67, 2.41, 3.14]
                ene_loss_max = [22. + 10, 28. + 10., 34. + 10., 40. + 10.]
                ene_loss_min = [22. - 10, 28. - 10., 34. - 10., 40. - 10.]
	if (theta_0_array[jexp] == 0.89):
		theta_exp = [1.68, 2.42, 3.15, 4.13]
                ene_loss_max = [24. + 10, 28. + 10., 32. + 10., 33. + 10.]
                ene_loss_min = [24. - 10, 28. - 10., 32. - 10., 33. - 10.]
	if (theta_0_array[jexp] == 1.06):
		theta_exp = [2.42, 3.15, 4.13]
                ene_loss_max = [38. + 10, 41. + 10., 43. + 10.]
                ene_loss_min = [38. - 10, 41. - 10., 43. - 10.]
	if (theta_0_array[jexp] == 1.23):
		theta_exp = [2.42, 3.15, 4.13]
                ene_loss_max = [38. + 10, 42. + 10., 44. + 10.]
                ene_loss_min = [38. - 10, 42. - 10., 44. - 10.]
    if (energy_0 == 500):
	if (theta_0_array[jexp] == 0.33):
		theta_exp = [0.66, 1.15, 1.64, 2.37, 3.11]
		ene_loss_max = [20. + 4, 32. + 4., 36. + 4., 38. + 4., 44. + 4.]
		ene_loss_min = [20. - 4, 32. - 4., 36. - 4., 38. - 4., 44. - 4.]
	if (theta_0_array[jexp] == 0.48):
		theta_exp = [1.15, 1.64, 2.37, 3.11]
		ene_loss_max = [32. + 4., 37. + 4., 38. + 4., 43. + 4.]
		ene_loss_min = [32. - 4., 37. - 4., 38. - 4., 43. - 4.]
	if (theta_0_array[jexp] == 0.64):
		theta_exp = [1.15, 1.64, 2.37, 3.11]
		ene_loss_max = [32. + 4., 35. + 4., 35. + 4., 38. + 4.]
		ene_loss_min = [32. - 4., 35. - 4., 35. - 4., 38. - 4.]
	if (theta_0_array[jexp] == 0.85):
		theta_exp = [1.65, 2.38, 3.12, 4.10]
		ene_loss_max = [34. + 4., 38. + 4., 37. + 4., 47. + 4.]
		ene_loss_min = [34. - 4., 38. - 4., 37. - 4., 47. - 4.]
	if (theta_0_array[jexp] == 1.02):
		theta_exp = [2.38, 3.12, 4.10]
		ene_loss_max = [35. + 4., 38. + 4., 42. + 4.]
		ene_loss_min = [35. - 4., 38. - 4., 42. - 4.]
	if (theta_0_array[jexp] == 1.19):
		theta_exp = [2.38, 3.12, 4.10]
		ene_loss_max = [35. + 4., 34. + 4., 40. + 4.]
		ene_loss_min = [35. - 4., 34. - 4., 40. - 4.]


    ax2.errorbar(theta_xrange_rebin, mean_energy_loss_rebin, xerr=err_theta_xrange_rebin, yerr=err_mean_energy_loss_rebin, capsize=0, fmt='o', lw = 2, color='black', ecolor=color_exp[jexp], label=str(theta_0_array[jexp])+' deg')
    ax2.fill_between(theta_exp, ene_loss_min, ene_loss_max, color='Gray', alpha=0.7, label='Diebold et al. (2015)')

    #ax2.set_yscale("log")
    ax2.set_xlabel("Scattering angle [deg.]")
    ax2.set_ylabel("Energy loss [keV]")
    if (energy_0 == 250): 
      ax2.set_ylim(0, 70)
    else:
      ax2.set_ylim(0, 300)
    ax2.set_xlim(0, 4.5)

    ax2.set_xlim(0, 4.5)

    title = str(N_in)+" events, "+str(energy_0)+" keV, "+str(theta_0_array[jexp])+" deg., "+model+" model"
    ax2.set_title(title)
    ax2.legend(numpoints=1, loc=2)

    plt.show()
    f_in.close()




