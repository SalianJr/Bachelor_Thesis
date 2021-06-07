# -*- coding: utf-8 -*-
"""
Created on Mon Mar 15 11:15:19 2021

@author: yashs
"""

# Imports
import pykep as pk
from pykep.orbit_plots import plot_planet, plot_lambert
from pykep import AU, DAY2SEC

# Plotting imports
import matplotlib as mpl
import matplotlib.pyplot as plt

# Lamberts Definition for E-J Trajectory
t1 = pk.epoch_from_string('2038-10-24 12:00:00') # Leaving Earth for Jupiter fly by
t2 = pk.epoch_from_string('2040-09-27 03:53:00') # Reaching Jupiter for Fly by
t4 = pk.epoch_from_string('2049-12-31 12:00:00') #Used for polotting purpose only, to keep the ephemerides well within the range

dt_EJ = (t2.mjd2000 - t1.mjd2000) * DAY2SEC

earth = pk.planet.jpl_lp('earth')
rE, vE = earth.eph(t1)

jupiter = pk.planet.jpl_lp('jupiter')
rJ, vJ = jupiter.eph(t2)

# Solving Lamberts problem for E-J Segment
l_EJ = pk.lambert_problem(r1 = rE, r2 = rJ, tof = dt_EJ, mu = pk.MU_SUN, max_revs=0)
v1_EJ = l_EJ.get_v1()[0]
v1_EJ = (v1_EJ[0]**2 + v1_EJ[1]**2 + v1_EJ[2]**2)**0.5   #magnitude of v1 while leaving earth- i.e start of trajectory
v2_EJ = l_EJ.get_v2()[0]
v2_EJ = (v2_EJ[0]**2 + v2_EJ[1]**2 + v2_EJ[2]**2)**0.5   #magnitude of v2 while arriving at jupiter- i.e end of trajectory   
DeltaV_EJ = abs(v2_EJ - v1_EJ)                        #delta V required to be on that trajectory


# Lamberts Definition for J-S Trajectory
t3 = pk.epoch_from_string('2045-03-17 12:00:00') # Reaching Saturn

dt_JS = (t3.mjd2000 - t2.mjd2000) * DAY2SEC

saturn = pk.planet.jpl_lp('saturn')
rS, vS = saturn.eph(t3)

# Solving Lamberts problem for J-S Segment
l_JS = pk.lambert_problem(r1 = rJ, r2 = rS, tof = dt_JS, mu = pk.MU_SUN, max_revs=0)
v1_JS = l_JS.get_v1()[0]
v1_JS = (v1_JS[0]**2 + v1_JS[1]**2 + v1_JS[2]**2)**0.5   #magnitude of v1 while leaving earth- i.e start of trajectory
v2_JS = l_JS.get_v2()[0]
v2_JS = (v2_JS[0]**2 + v2_JS[1]**2 + v2_JS[2]**2)**0.5   #magnitude of v2 while arriving at saturn- i.e end of trajectory   
DeltaV_JS = abs(v2_JS - v1_JS)                        #delta V required to be on that trajectory

Total_DV = DeltaV_JS + DeltaV_EJ


# We plot
mpl.rcParams['legend.fontsize'] = 11

# Create the figure and axis for E-J
fig1 = plt.figure(figsize = (23,10))

ax1 = fig1.add_subplot(1, 3, 2, projection='3d')
ax1.scatter([0], [0], [0], color=['y'], s=50)
ax1.view_init(90, 90)

ax2 = fig1.add_subplot(1, 3, 3, projection='3d')
ax2.scatter([0], [0], [0], color=['y'], s=50)
ax2.view_init(90,90)
    
# Axes 1 planet positions
plot_planet(earth, t0=t1, color=(0, 0.7, .9), legend=True, units=AU, axes=ax1)
plot_planet(jupiter, t0=t1, tf=t4, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax1)

# Axes 2 with lamberts Trajectory
plot_planet(earth, t0=t1, color=(0, 0.7, .9), legend=True, units=AU, axes=ax2)
plot_planet(jupiter, t0=t2, tf=t4, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax2)
axis = plot_lambert(l_EJ, color='b', legend=True, units=AU, axes=ax2) 

# Create the figure and axis for J-S
fig2 = plt.figure(figsize = (23,10))

ax1 = fig2.add_subplot(1, 3, 2, projection='3d')
ax1.scatter([0], [0], [0], color=['y'], s=50)
ax1.view_init(90, 90)

ax2 = fig2.add_subplot(1, 3, 3, projection='3d')
ax2.scatter([0], [0], [0], color=['y'], s=50)
ax2.view_init(90,90)

# Axes 1 planet positions
plot_planet(jupiter, t0=t2, tf=t4, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax1)
plot_planet(saturn, t0=t2, tf=t4, color=(.67, .37, 0.29), s=150, legend=True, units=AU, axes=ax1)

# Axes 2 with lamberts Trajectory
plot_planet(jupiter, t0=t2, tf=t4, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax2)
plot_planet(saturn, t0=t3, tf=t4, color=(.67, .37, 0.29), s=150, legend=True, units=AU, axes=ax2)
axis = plot_lambert(l_JS, color='b', legend=True, units=AU, axes=ax2) 

# Create the figure and axis for E-J-S
fig3 = plt.figure(figsize = (30,14))

ax1 = fig3.add_subplot(1, 3, 2, projection='3d')
ax1.scatter([0], [0], [0], color=['y'], s=50)
ax1.view_init(90, 90)

# E-J
plot_planet(earth, t0=t1, color=(0, 0.7, .9), legend=True, units=AU, axes=ax1)
plot_planet(jupiter, t0=t2, tf=t4, color=(1, .6, 0), s=300, legend=True, units=AU, axes=ax1)
axis = plot_lambert(l_EJ, color='b', legend=True, units=AU, axes=ax1)

# J-S
plot_planet(jupiter, t0=t2, tf=t4, color=(1, .6, 0), s=300, units=AU, axes=ax1)
plot_planet(saturn, t0=t3, tf=t4, color=(.67, .37, 0.29), s=150, legend=True, units=AU, axes=ax1)
axis = plot_lambert(l_JS, color='b', units=AU, axes=ax1)
