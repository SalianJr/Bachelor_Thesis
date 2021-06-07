# -*- coding: utf-8 -*-
"""
Created on Tue Dec 15 10:29:46 2020

@author: yashs
"""
# Imports
import pykep as pk
from pykep.orbit_plots import plot_planet, plot_lambert
from pykep import AU, DAY2SEC

# Plotting imports
import matplotlib as mpl
import matplotlib.pyplot as plt

# We define the Lambert problem
t1 = pk.epoch_from_string('2011-08-05 16:25:00')
t2 = pk.epoch_from_string('2016-07-05 03:53:00')
dt = (t2.mjd2000 - t1.mjd2000) * DAY2SEC

earth = pk.planet.jpl_lp('earth')
rE, vE = earth.eph(t1)

jupiter = pk.planet.jpl_lp('jupiter')
rJ, vJ = jupiter.eph(t2)

# We solve the Lambert problem
l = pk.lambert_problem(r1 = rE, r2 = rJ, tof = dt, mu = pk.MU_SUN, max_revs=0)
v1 = l.get_v1()[0]
v1 = (v1[0]**2 + v1[1]**2 + v1[2]**2)**0.5   #magnitude of v1 while leaving earth- i.e start of trajectory
v2 = l.get_v2()[0]
v2 = (v2[0]**2 + v2[1]**2 + v2[2]**2)**0.5   #magnitude of v2 while arriving at jupiter- i.e end of trajectory   
DeltaV = abs(v2 - v1)                        #delta V required to be on that trajectory

# We plot
mpl.rcParams['legend.fontsize'] = 10

# Create the figure and axis
fig = plt.figure(figsize = (20,7))

ax1 = fig.add_subplot(1, 3, 2, projection='3d')
ax1.scatter([0], [0], [0], color=['y'])
ax1.view_init(90, 90)

ax2 = fig.add_subplot(1, 3, 3, projection='3d')
ax2.scatter([0], [0], [0], color=['y'])
ax2.view_init(90,90)

# for ax in [ax1, ax2]:
#     # Plot the planet orbits
#     plot_planet(earth, t0=t1, color=(0, 0.7, .9), legend=True, units=AU, axes=ax)
#     plot_planet(jupiter, t0=t2, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax)
  
    # # Plot the Lambert solutions
    # axis = plot_lambert(l, color='b', legend=True, units=AU, axes=ax)
    
plot_planet(earth, t0=t1, color=(0, 0.7, .9), legend=True, units=AU, axes=ax1)
plot_planet(jupiter, t0=t1, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax1)
#axis = plot_lambert(l, color='b', legend=True, units=AU, axes=ax1)

plot_planet(earth, t0=t1, color=(0, 0.7, .9), legend=True, units=AU, axes=ax2)
plot_planet(jupiter, t0=t2, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax2)
axis = plot_lambert(l, color='b', legend=True, units=AU, axes=ax2)      
    
    
    
    
    
    