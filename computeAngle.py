import numpy as np
theta_0 = 1.17
print "theta_0 = ", theta_0

theta_0_90 = 90. + theta_0   # deg.
theta_0_rad = theta_0*(np.pi/180.)
phi_0 = 0.

Px = - np.sin(theta_0_rad)*np.cos(phi_0)
Py = - np.sin(theta_0_rad)*np.sin(phi_0)
Pz = - np.cos(theta_0_rad)

print "Px, Py, Pz = ", Px, Py, Pz


height = 5.5*np.tan(theta_0*(np.pi/180.))
print "height [cm]: ", height

scatt_angle = 3.08 - theta_0 #deg.
d = 933. #mm
print "-------------------------"
print "Scattering angle [deg]: ", scatt_angle
print "radial distance [mm]: ", d
print "Z center: ", np.sin(scatt_angle*(np.pi/180.))*d
print "X distance: ", np.cos(scatt_angle*(np.pi/180.))*d


