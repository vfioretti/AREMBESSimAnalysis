"""
 plotTotalScatteringEfficiency_50_1000.py  -  description
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
rebin = int(arg_list[2])

# experimental data
# 1MeV
N_exp = 6
color_exp = ['red', 'green', 'blue', 'purple', 'orange']
theta_0_array = [0.3, 0.46, 0.61, 0.89, 1.0, 1.17]

# Plot the results
fig = plt.figure(1, figsize=(15, 10))
ax = fig.add_subplot(111)

handles2 = []

models = ['SingleScattering', 'SingleScattering']
models_title = ['E$_{0}$ = 50 keV', 'E$_{0}$ = 1000 keV']
energy_0_array = [50, 1000]

model_fmt = ['-o', '-^']

for jexp in xrange(N_exp):
 for jm in xrange(len(models)):
    
    # Plot the results
    fig2 = plt.figure(1, figsize=(15, 10))
    ax2 = fig2.add_subplot(111)

    # Reading simulation output
    name_file_in = './diebold_test_results/'+str(int(energy_0_array[jm]))+'keV_'+str(theta_0_array[jexp])+'_'+str(N_in)+'/'+models[jm]+'.dat'
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


 ax2.set_yscale("log")
 ax2.set_xlabel("Scattering angle [deg.]")
 ax2.set_ylabel("Scattering efficiency [sr$^{-1}$]")
 ax2.set_ylim(10, 100000)
 ax2.set_xlim(0, 4.5)
 title = " Single scattering, incident angle = "+str(theta_0_array[jexp])+" deg."
 ax2.set_title(title)
 ax2.legend(numpoints=1, loc=1)

 plt.grid()
 plt.show()

