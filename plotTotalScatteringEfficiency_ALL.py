"""
 plotTotalScatteringEfficiency_ALL.py  -  description
 ---------------------------------------------------------------------------------
 Plotting the soft proton scattering efficiency
 ---------------------------------------------------------------------------------
 copyright            : (C) 2016 Valentina Fioretti
 email                : fioretti@iasfbo.inaf.it
 ----------------------------------------------
 Usage:
 python plotTotalScatteringEfficiency_ALL.py N_in energy_0 rebin
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

# experimental data
# 1MeV
N_exp = 6
color_exp = ['red', 'green', 'blue', 'purple', 'orange']
if energy_0 == 1000:
    theta_0_array = [0.3, 0.46, 0.61, 0.89, 1.0, 1.17]
if energy_0 == 500:
    theta_0_array = [0.33, 0.48, 0.64, 0.85, 1.02, 1.19]
if energy_0 == 250:
    theta_0_array = [0.36, 0.51, 0.67, 0.89, 1.06, 1.23]

# Plot the results
fig = plt.figure(1, figsize=(15, 10))
ax = fig.add_subplot(111)

handles2 = []

#models = ['Firsov', 'Remizovich', 'MultipleScattering']
#models_title = ['Firsov', 'Remizovich (elastic)', 'Multiple Scattering']
models = ['Firsov', 'Remizovich', 'MultipleScattering', 'MultipleScattering_Option4', 'SingleScattering']
models_title = ['Firsov', 'Remizovich (elastic)', 'Multiple Scattering (option 3)', 'Multiple Scattering (option 4)', 'Single Scattering']


model_fmt = ['-o', '-^', '-s', '-<', '->']

for jexp in xrange(N_exp):
 for jm in xrange(len(models)):
    
    # Plot the results
    fig2 = plt.figure(1, figsize=(15, 10))
    ax2 = fig2.add_subplot(111)

    # Reading simulation output
    name_file_in = './diebold_test_results/'+str(int(energy_0))+'keV_'+str(theta_0_array[jexp])+'_'+str(N_in)+'/'+models[jm]+'.dat'
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
                err_angle_tot+=err_angle[jel + i]

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


    ax2.errorbar(angle_array_rebin, N_array_out_norm_rebin, xerr=err_angle_rebin, yerr=err_N_array_out_rebin, capsize=0, fmt=model_fmt[jm], ms=10, lw = 2, ecolor='black', color=color_exp[jm], label='BoGEMMS - '+models_title[jm])

    f_in.close()


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
 if (energy_0 == 250):
	if (theta_0_array[jexp] == 0.36):
		theta_exp = [0.69, 1.18, 1.67, 2.41, 3.14]
		eff_exp_max = [340.11 + 88.15, 468.94 + 121.65, 224.33 + 58.22, 84.03 + 21.83, 42.45 + 11.07]
		eff_exp_min = [340.11 - 88.15, 468.94 - 121.65, 224.33 - 58.22, 84.03 - 21.83, 42.45 - 11.07]
	if (theta_0_array[jexp] == 0.51):
		theta_exp = [0.69, 1.18, 1.67, 2.41, 3.14]
		eff_exp_max = [1.24 + 0.52, 311.78 + 80.94, 231.58 + 60.12, 99.62 + 25.87, 51.51 + 13.45]
		eff_exp_min = [1.24 - 0.52, 311.78 - 80.94, 231.58 - 60.12, 99.62 - 25.87, 51.51 - 13.45]
	if (theta_0_array[jexp] == 0.67):
		theta_exp = [1.18, 1.67, 2.41, 3.14]
		eff_exp_max = [139.60 + 36.28, 189.61 + 49.23, 107.04 + 27.81, 56.75 + 14.79]
		eff_exp_min = [139.60 - 36.28, 189.61 - 49.23, 107.04 - 27.81, 56.75 - 14.79]
	if (theta_0_array[jexp] == 0.89):
		theta_exp = [1.68, 2.42, 3.15, 4.13]
		eff_exp_max = [161.91 + 42.01, 114.94 + 29.82, 64.40 + 16.72, 30.44 + 7.93]
		eff_exp_min = [161.91 - 42.01, 114.94 - 29.82, 64.40 - 16.72, 30.44 - 7.93]
	if (theta_0_array[jexp] == 1.06):
		theta_exp = [2.42, 3.15, 4.13]
		eff_exp_max = [108.04 + 28.06, 66.12 + 17.17, 32.11 + 8.36]
		eff_exp_min = [108.04 - 28.06, 66.12 - 17.17, 32.11 - 8.36]
	if (theta_0_array[jexp] == 1.23):
		theta_exp = [2.42, 3.15, 4.13]
		eff_exp_max = [96.72 + 25.13, 69.99 + 18.16, 36.33 + 9.44]
		eff_exp_min = [96.72 - 25.13, 69.99 - 18.16, 36.33 - 9.44]
 if (energy_0 == 500):
	if (theta_0_array[jexp] == 0.33):
		theta_exp = [0.66, 1.15, 1.64, 2.37, 3.11]
		eff_exp_max = [436.13 + 113.27, 578.17 + 149.86, 252.38 + 65.51, 94.10 + 24.44, 45.47 + 11.82]
		eff_exp_min = [436.13 - 113.27, 578.17 - 149.86, 252.38 - 65.51, 94.10 - 24.44, 45.47 - 11.82]
	if (theta_0_array[jexp] == 0.48):
		theta_exp = [0.66, 1.15, 1.64, 2.37, 3.11]
		eff_exp_max = [0.20 + 0.57, 388.00 + 100.53, 262.58 + 68.16, 108.76 + 28.26, 55.08 + 14.32]
		eff_exp_min = [0.20 - 0.57, 388.00 - 100.53, 262.58 - 68.16, 108.76 - 28.26, 55.08 - 14.32]
	if (theta_0_array[jexp] == 0.64):
		theta_exp = [1.15, 1.64, 2.37, 3.11]
		eff_exp_max = [142.36 + 36.98, 226.93 + 58.93, 130.30 + 33.86, 68.80 + 17.87]
		eff_exp_min = [142.36 - 36.98, 226.93 - 58.93, 130.30 - 33.86, 68.80 - 17.87]
	if (theta_0_array[jexp] == 0.85):
		theta_exp = [1.65, 2.38, 3.12, 4.10]
		eff_exp_max = [181.05 + 47.01, 124.62 + 32.36, 66.06 + 17.18, 32.23 + 8.46]
		eff_exp_min = [181.05 - 47.01, 124.62 - 32.36, 66.06 - 17.18, 32.23 - 8.46]
	if (theta_0_array[jexp] == 1.02):
		theta_exp = [2.38, 3.12, 4.10]
		eff_exp_max = [116.66 + 30.30, 68.31 + 17.77, 36.72 + 9.62]
		eff_exp_min = [116.66 - 30.30, 68.31 - 17.77, 36.72 - 9.62]
	if (theta_0_array[jexp] == 1.19):
		theta_exp = [2.38, 3.12, 4.10]
		eff_exp_max = [103.50 + 26.88, 70.69 + 18.39, 36.90 + 9.67]
		eff_exp_min = [103.50 - 26.88, 70.69 - 18.39, 36.90 - 9.67]

 ax2.fill_between(theta_exp, eff_exp_min, eff_exp_max, color='Gray', alpha=0.7, label='Diebold et al. (2015)')

 ax2.set_yscale("log")
 ax2.set_xlabel("Scattering angle [deg.]")
 ax2.set_ylabel("Scattering efficiency [sr$^{-1}$]")
 ax2.set_ylim(10, 100000)
 ax2.set_xlim(0, 4.5)
 title = " E$_{0}$ = "+str(energy_0)+" keV, incident angle = "+str(theta_0_array[jexp])+" deg."
 ax2.set_title(title)
 ax2.legend(numpoints=1, loc=1)

 plt.grid()
 plt.show()

#ax.set_yscale("log")
#ax.set_xlabel("Scattering angle [deg.]")
#ax.set_ylabel("Scattering efficiency [sr$^{-1}$]")

#title = str(N_in)+" events, "+str(energy_0)+" keV, "+model+" model"
#ax.set_title(title)

#first_legend = plt.legend(handles=[legend11, legend12], numpoints=1)
#ax = plt.gca().add_artist(first_legend)
#second_legend = plt.legend(handles=handles2, numpoints=1, loc=4)


