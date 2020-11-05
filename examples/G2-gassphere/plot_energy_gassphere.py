#!/usr/bin/env python3
"""
Code that plots the energy evolution for the Evrard collapse

"""
# load libraries
import numpy as np
import matplotlib.pyplot as plt
import csv


time = np.zeros(61)
thermal_energy = np.zeros(61)
potential_energy = np.zeros(61)
kinetic_energy = np.zeros(61)

with open('output/energy.txt') as csvfile:
    readCSV = csv.reader(csvfile)
    line = 0
    for row in readCSV:
        values = row[0].split()
        if len(values) != 8:
            break
        time[line] = values[0]
        thermal_energy[line] = values[1]
        potential_energy[line] = values[2]
        kinetic_energy[line] = values[3]


        line = line + 1

total_energy = thermal_energy + potential_energy + kinetic_energy

plt.plot(time, thermal_energy,label="Thermal energy")
plt.plot(time, potential_energy,label="Potential energy")
plt.plot(time, kinetic_energy,label="Kinetic energy")
plt.plot(time, total_energy,label="Total energy")
plt.plot(time, (total_energy-total_energy[0])/total_energy[0])
plt.ylabel('Energy')
plt.xlabel('Time')
plt.legend()
plt.show()