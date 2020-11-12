#!/usr/bin/env python3
"""
Code that plots different radial profiles for the Evrard collapse.
The results are compared with 1D PPM results from (Steinmetz & Mueller 1993) for t= 0.8

"""
# load libraries
import sys  # load sys; needed for exit codes
import numpy as np  # load numpy
import h5py  # load h5py; needed to read snapshots
import matplotlib
import matplotlib.pyplot as plt  ## needs to be active for plotting!
import csv

matplotlib.rc_file_defaults()
FloatType = np.float64

#loads the gadget snapshot name
def load_snapshot(name):
    gamma = 5./3.

    try:
        data = h5py.File(name, "r")
    except:
        print("could not open file : "+ name+" !")
        exit(1)

    time = FloatType(data["Header"].attrs["Time"])
    Pos = np.array(data["PartType0"]["Coordinates"], dtype = FloatType) 
    Density = np.array(data["PartType0"]["Density"], dtype = FloatType)
    Velocity = np.array(data["PartType0"]["Velocities"], dtype = FloatType)
    Uthermal = np.array(data["PartType0"]["InternalEnergy"], dtype = FloatType)

    radius = np.sqrt(Pos[:,0]**2 + Pos[:,1]**2 + Pos[:,2]**2)
    vr = (Velocity[:,0] * Pos[:,0] + Velocity[:,1] * Pos[:,1] + Velocity[:,2] * Pos[:,2]) / radius[:]
    origin_particles = np.argwhere(radius == 0.0)
    vr[origin_particles] = 0

    print(np.max(Density))
    press = (gamma-1)*Uthermal*Density

    A = press / Density**gamma

    return radius,Density, vr, A, time 

#bins the particle properties
def get_data_bins(radius,Density, vr, A):
	number_bins = 100
	min_r = 0.01
	max_r = 1.0

	r_bin = np.zeros(number_bins)
	vr_bin = np.zeros(number_bins)

	rho_bin = np.zeros(number_bins)
	count_bin = np.zeros(number_bins)
	entropy_bin = np.zeros(number_bins)

	for i in range(radius.size):
		if(radius[i] < min_r or radius[i] > max_r):
			continue
		bin_value = int((np.log10(radius[i] / min_r))/(np.log10(max_r/min_r)) * number_bins)
		count_bin[bin_value] = count_bin[bin_value] + 1
		vr_bin[bin_value] = vr_bin[bin_value] + vr[i]
		rho_bin[bin_value] = rho_bin[bin_value] + Density[i]
		entropy_bin[bin_value] = entropy_bin[bin_value] + A[i]


	vr_bin /= count_bin
	rho_bin /= count_bin
	entropy_bin /= count_bin

	for i in range(number_bins):
		r_bin[i] = (i+0.5)* (np.log10(max_r/min_r)/number_bins) + np.log10(min_r)
		r_bin[i] = 10**r_bin[i]
		print(count_bin[i])

	return r_bin,rho_bin, vr_bin, entropy_bin

#loads 1D PPM results
def load_ppm_result():
	gamma = 5./3.
	rost = 3./4./np.pi
	est = 1.054811e-1  / 1.05
	pst = rost*est
	vst = np.sqrt(est)
	rst = 1.257607
	time = 0

	radius = np.zeros(350)
	rho = np.zeros(350)
	vr = np.zeros(350)
	press = np.zeros(350)

	with open('ppm_profile/ppm1oaf') as csvfile:
		readCSV = csv.reader(csvfile)
		line = 0
		for row in readCSV:
			line = line+1
			values = row[0].split()
			if(line == 1):
				time = values[1]
				continue
			if(line == 352):
				break

			radius[line -2] = float(values[1]) /rst*1e-11
			rho[line -2] = float(values[2]) /rost
			vr[line -2] = float(values[4]) /vst*1e-8
			press[line -2] = float(values[3])/pst*1e-16

	rho=rho*(3.0/(4*np.pi))
	press = press*(3.0/(4*np.pi))

	entropy = press / rho**gamma

	return radius, rho, vr, entropy


i_file = int(sys.argv[1])

filename = 'output/snapshot_%03d.hdf5' % i_file
    
radius,Density, vr, A, time = load_snapshot(filename)
    
radius,Density, vr, A = get_data_bins(radius,Density, vr, A)

fig, ax = plt.subplots(1,3,figsize=(18,6))

    
ax[0].scatter(radius, Density, facecolors='none', edgecolors='b',alpha = 0.5)
ax[1].scatter(radius, vr, facecolors='none', edgecolors='b',alpha = 0.5)
ax[2].scatter(radius, A, facecolors='none', edgecolors='b',alpha = 0.5)

radius_ppm, rho_ppm, vr_ppm, entropy_ppm = load_ppm_result()

ax[0].plot(radius_ppm, rho_ppm, color="k")
ax[1].plot(radius_ppm, vr_ppm, color="k")
ax[2].plot(radius_ppm, entropy_ppm, color="k")

ax[0].set_xscale('log')
ax[1].set_xscale('log')
ax[2].set_xscale('log')
ax[0].set_yscale('log')
ax[0].set_xlim(0.01,1)
ax[1].set_xlim(0.01,1)
ax[2].set_xlim(0.01,1)

ax[0].set_ylim(0.004,700)
ax[1].set_ylim(-1.8,0)
ax[2].set_ylim(0.0,0.2)

ax[0].set_xlabel("R")
ax[1].set_xlabel("R")
ax[2].set_xlabel("R")

ax[0].set_ylabel(r"$\rho$")
ax[1].set_ylabel(r"$V_R$")
ax[2].set_ylabel(r"$P/\rho^\gamma$")
plt.tight_layout()




plt.savefig("evrard_"+str(i_file)+".eps")
plt.show()

