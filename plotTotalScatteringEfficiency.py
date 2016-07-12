"""
 plotTotalScatteringEfficency.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the soft proton scattering efficiency
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotTotalScatteringEfficency.py filedir N_in energy_0 model
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
model = arg_list[3]

# experimental data
# 1MeV
N_exp = 6
rebin = 5
color_exp = ['Gray', 'red', 'green', 'blue', 'purple', 'OrangeRed']
theta_0_array = [0.3, 0.46, 0.61, 0.89, 1.0, 1.17]

# Plot the results
fig = plt.figure(1, figsize=(15, 10))
ax = fig.add_subplot(111)

handles2 = []

for jexp in xrange(N_exp):
    
    # Plot the results
    fig2 = plt.figure(1, figsize=(15, 10))
    ax2 = fig2.add_subplot(111)

    # Reading simulation output
    name_file_in = './diebold_test_results/'+str(int(energy_0))+'keV_'+str(theta_0_array[jexp])+'_'+str(N_in)+'/'+model+'.dat'
    f_in = open(name_file_in, 'r')
	
    angle_array = []
    left_angle_array = []
    err_angle = []
    N_array_out_norm = []
    err_N_array_out = []
    N_array_out = []
    N_in_eff = []
    solid_angle_array = []

    angle_array_rebin = []
    err_angle_rebin = []
    N_array_out_norm_rebin = []
    err_N_array_out_rebin = []
    N_array_out_rebin = []
    N_in_eff_rebin = []
    solid_angle_array_rebin = []

    for line in f_in:
        line = line.strip()
        columns = line.split()
        # angle_x err_angle_x Eff err_Eff N_out N_in solid_angle
        columns[0] = float(columns[0])
        columns[1] = float(columns[1])
        columns[2] = float(columns[2])
        columns[3] = float(columns[3])
        columns[4] = float(columns[4])
        columns[5] = float(columns[5])
        columns[6] = float(columns[6])
        
        angle_array.append(columns[0])
        err_angle.append(columns[1])
        N_array_out_norm.append(columns[2])
        err_N_array_out.append(columns[3])
        N_array_out.append(columns[4])
        N_in_eff.append(columns[5])
        solid_angle_array.append(columns[6])


    jel=0
    while(1):
        if (jel + rebin) < len(angle_array):
            N_out_tot = 0
            N_in_tot = 0
            solid_angle_tot = 0
            err_angle_tot = 0
            for i in xrange(rebin):
                N_out_tot+=N_array_out[jel + i]
                N_in_tot=N_in_eff[jel + i]
                solid_angle_tot+=solid_angle_array[jel + i]
                err_angle_tot+=err_angle[i]

            angle_array_rebin.append(angle_array[jel] + (angle_array[jel+ rebin-1] - angle_array[jel])/2.)
            err_angle_rebin.append(err_angle_tot)
            N_array_out_rebin.append(N_out_tot)
            N_in_eff_rebin.append(N_in_tot)
            solid_angle_array_rebin.append(solid_angle_tot)
        else:
            break
        jel = jel + rebin

    for jel in xrange(len(angle_array_rebin)):
         N_array_out_norm_rebin.append((N_array_out_rebin[jel]/N_in_eff_rebin[jel])/solid_angle_array_rebin[jel])
         err_N_array_out_rebin.append((np.sqrt(N_array_out_rebin[jel])/N_in_eff_rebin[jel])/solid_angle_array_rebin[jel])


    if (energy_0 == 1000):
        if (theta_0_array[jexp] == 0.3):
            theta_exp = [0.63, 1.12, 1.61, 2.35, 3.08]
            eff_exp_max = [415.06 + 107.69, 584.64 + 151.61, 250.47 + 65.02, 98.88 + 25.71, 41.60 + 10.84]
            eff_exp_min = [415.06 - 107.69, 584.64 - 151.61, 250.47 - 65.02, 98.88 - 25.71, 41.60 - 10.84]
        if (theta_0_array[jexp] == 0.46):
            theta_exp = [0.63, 1.12, 1.61, 2.35, 3.08]
            eff_exp_max = [2.09 + 1.25, 413.98 + 107.47, 301.39 + 78.35, 116.20 + 30.21, 57.49 + 14.96]
            eff_exp_min = [2.09 - 1.25, 413.98 - 107.47, 301.39 - 78.35, 116.20 - 30.21, 57.49 - 14.96]
        if (theta_0_array[jexp] == 0.61):
            theta_exp = [1.12, 1.61, 2.35, 3.08]
            eff_exp_max = [160.94 + 41.89, 244.21 + 63.41, 136.86 + 35.61, 65.08 + 16.93]
            eff_exp_min = [160.94 - 41.89, 244.21 - 63.41, 136.86 - 35.61, 65.08 - 16.93]
        if (theta_0_array[jexp] == 0.89):
            theta_exp = [1.68, 2.42, 3.15, 4.13]
            eff_exp_max = [161.91 + 42.01, 114.94 + 29.82, 64.40 + 16.72, 30.44 + 7.93]
            eff_exp_min = [161.91 - 42.01, 114.94 - 29.82, 64.40 - 16.72, 30.44 - 7.93]
        if (theta_0_array[jexp] == 1.00):
            theta_exp = [2.35, 3.09, 4.07]
            eff_exp_max = [132.88 + 34.53, 80.00 + 20.80, 37.47 + 9.76]
            eff_exp_min = [132.88 - 34.53, 80.00 - 20.80, 37.47 - 9.76]
        if (theta_0_array[jexp] == 1.17):
            theta_exp = [2.35, 3.09, 4.07]
            eff_exp_max = [113.94 + 29.62, 84.25 + 21.91, 38.90 + 10.13]
            eff_exp_min = [113.94 - 29.62, 84.25 - 21.91, 38.90 - 10.13]

    """
    if jexp == 0:
       legend11 = ax.errorbar(angle_array_rebin, N_array_out_norm_rebin, xerr=err_angle_rebin, yerr=err_N_array_out_rebin, capsize=0, fmt='-o', lw = 2, color='black', ecolor='black', label='BoGEMMS')
       legend12 = ax.fill_between(theta_exp, eff_exp_min, eff_exp_max, color='black', label='Experiment')
       handles2.append(ax.errorbar(angle_array_rebin, N_array_out_norm_rebin, xerr=err_angle_rebin, yerr=err_N_array_out_rebin, capsize=0, fmt='-o', lw = 2, color=color_exp[jexp], ecolor=color_exp[jexp], label = str(theta_0_array[jexp])+' deg.'))
       ax.fill_between(theta_exp, eff_exp_min, eff_exp_max, color=color_exp[jexp], alpha=0.7)

    else:
       handles2.append(ax.errorbar(angle_array_rebin, N_array_out_norm_rebin, xerr=err_angle_rebin, yerr=err_N_array_out_rebin, capsize=0, fmt='-o', lw = 2, color=color_exp[jexp], ecolor=color_exp[jexp], label = str(theta_0_array[jexp])+' deg.'))
       ax.fill_between(theta_exp, eff_exp_min, eff_exp_max, color=color_exp[jexp], alpha=0.7)
    """

    ax2.errorbar(angle_array_rebin, N_array_out_norm_rebin, xerr=err_angle_rebin, yerr=err_N_array_out_rebin, capsize=0, fmt='o', lw = 2, color='black', ecolor=color_exp[jexp], label='BoGEMMS')
    ax2.fill_between(theta_exp, eff_exp_min, eff_exp_max, color=color_exp[jexp], alpha=0.7, label='Experiment')

    ax2.set_yscale("log")
    ax2.set_xlabel("Scattering angle [deg.]")
    ax2.set_ylabel("Scattering efficiency [sr$^{-1}$]")
    ax2.set_ylim(1, 10000)
    ax2.set_xlim(0, 4.5)
    title = str(N_in)+" events, "+str(energy_0)+" keV, "+str(theta_0_array[jexp])+" deg., "+model+" model"
    ax2.set_title(title)
    ax2.legend(numpoints=1, loc=1)

    plt.show()
    f_in.close()


#ax.set_yscale("log")
#ax.set_xlabel("Scattering angle [deg.]")
#ax.set_ylabel("Scattering efficiency [sr$^{-1}$]")

#title = str(N_in)+" events, "+str(energy_0)+" keV, "+model+" model"
#ax.set_title(title)

#first_legend = plt.legend(handles=[legend11, legend12], numpoints=1)
#ax = plt.gca().add_artist(first_legend)
#second_legend = plt.legend(handles=handles2, numpoints=1, loc=4)


