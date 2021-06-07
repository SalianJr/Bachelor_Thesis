# -*- coding: utf-8 -*-
"""
Version 1.4
Created on Sat Jun  5 12:24:04 2021

@author: yashs
"""

# Fly-by trajectories using Pykep and Pygmo
# UDP created using mga_1dsm class of pykep
# Population with iteration, Island, Archipelago - all variations tried out
# Multi-objective Implemented on two islands for the given udp
# Made the code more compact by creating run_island function, to instantiate multiple islands and present the optimized results in the form of graphs and text
# Added another functionality to plot all the evolved populations of instantiated islands using scatter plots custom function

import pykep as pk
import pygmo as pg
import numpy as np

# Timing Imports
from timeit import default_timer as timer
from datetime import timedelta

# Plotting Imports
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
import matplotlib.pyplot as plt

# Create, Evolve and plot multiple Islands together.
def run_isl(isl,archi,udp,evl=20):   #returns tuple with number of (single,multi) objective islands
    j = 1 
    single = 0
    multi= 0
    for i in isl:
        if i.get_population().get_f()[0].size == 1:     # Single objective island
            print()    
            print("------------------------------------------------Island ",j,"-----------------------------------------------")
            print(i.get_algorithm())
            print("-------------------------------------------------------------------------------------------------------")

            start = timer()

            i.evolve(evl)
            i.wait()

            end = timer()

            f_sol = i.get_population().champion_f
            x_sol = i.get_population().champion_x

            print("############################### Champion Function Values ###########################################")
            print(f_sol)
            print()
            print("############################### Champion Decision Vectors ###########################################")
            print(x_sol)
            print()
            print('Time for implementation: ',timedelta(seconds=end-start))
            print()
            archi[0].push_back(i)
            udp.pretty(x_sol)
            udp.plot(x_sol)
            j = j+1
            single = single+1
            
        else:  #Multi-Objective 
            print()
            print("----------------------------------------------Island ",j,"-----------------------------------------------")
            print(i.get_algorithm())
            print("-------------------------------------------------------------------------------------------------------")

            start = timer()

            i.evolve(evl)
            i.wait()

            end = timer()

            sols_f = i.get_population().get_f()
            sols_x = i.get_population().get_x()
            f_sol_T = sols_f[np.argmin(np.hsplit(sols_f,2)[1])]
            x_sol_T = sols_x[np.argmin(np.hsplit(sols_f,2)[1])]
            f_sol_V = sols_f[np.argmin(np.hsplit(sols_f,2)[0])]
            x_sol_V = sols_x[np.argmin(np.hsplit(sols_f,2)[0])]

            print("############################### Best Time Value and Vector ###########################################")
            print(f_sol_T,'\n')
            print(x_sol_T)
            print()
            print("############################### Best Delta V Value and Vector ###########################################")
            print(f_sol_V,'\n')
            print(x_sol_V)
            print()
            print('Time for implementation: ',timedelta(seconds=end-start))
            print()
            archi[1].push_back(i)
            j = j+1
            multi = multi+1
    return (single,multi)

# Evolving Single and Multi-Objective Archipelagos - isl_type_num is returned from run_isl which can be passed here
def run_archi(archi,udp,isl_type_num =(1,1),evl=20):
    j = 1
    for a in archi:
        if a[0].get_population().get_f()[0].size == 1:
            print(a)
            start = timer()

            a.evolve(evl)
            a.wait()

            end = timer()

            f_sols = a.get_champions_f()
            x_sols = a.get_champions_x()
            print("----------------------------------------------Archipelago ",j,"-----------------------------------------------")
            print("############################### Champion Function Values ###########################################")
            print(f_sols)
            print()
            print("############################### Champion Decision Vectors ###########################################")
            print(x_sols)
            print()
            print('Time for implementation: ',timedelta(seconds=end-start))
            print()
            print("############################### Trajectory Info ###########################################")
            idx = f_sols.index(min(f_sols))
            udp.pretty(a.get_champions_x()[idx])
            udp.plot(a.get_champions_x()[idx])
            j = j+1

        else:
            print(a)

            start = timer()

            a.evolve(evl)
            a.wait()

            end = timer()
            
            sols_f_T = list()
            sols_f_V = list()
            idx_min_T = list()
            idx_min_V = list()
            
            for i in range(isl_type_num[1]):
                archi_mo_sols_f = a[i].get_population().get_f()
                idx_min_T = idx_min_T + [np.argmin(np.hsplit(archi_mo_sols_f,2)[1])]
                sols_f_T = sols_f_T + [np.hsplit(archi_mo_sols_f,2)[1][idx_min_T[i]]]
                idx_min_V = idx_min_V + [np.argmin(np.hsplit(archi_mo_sols_f,2)[0])]
                sols_f_V = sols_f_V + [np.hsplit(archi_mo_sols_f,2)[0][idx_min_V[i]]]
                
            champion_x_T = a[np.argmin(np.array(sols_f_T))].get_population().get_x()[idx_min_T[np.argmin(np.array(sols_f_T))]]
            champion_x_V = a[np.argmin(np.array(sols_f_V))].get_population().get_x()[idx_min_V[np.argmin(np.array(sols_f_V))]]
            
            print('Time for implementation: ',timedelta(seconds=end-start))
            print()
            print("############################### Best Time Trajectory Info ###########################################")
            print(min(sols_f_T))
            print()
            print(champion_x_T)
            print()
            udp.pretty(champion_x_T)
            udp.plot(champion_x_T)
            print()
            
            print("############################### Best Velocity Trajectory Info ###########################################")
            print(min(sols_f_V))
            print()
            print(champion_x_V)
            print()
            udp.pretty(champion_x_V)
            udp.plot(champion_x_V)
            print()  
            
# Scatter plot function for all single and multi-objective optimizations
def scatter_all_sols(isl, add_vinf_arr=False, marker_size=10, dpi=500, fig_size=(10.87402,7.8663)):
    j = 1
    for i in isl:
        if i.get_population().get_f()[0].size == 2:
    
            x = np.hsplit(i.get_population().get_f(),2)[1] / 365.25    #years
            y = np.hsplit(i.get_population().get_f(),2)[0] / 1000  #km/s delta v
            
            if add_vinf_arr==True:
                _ = np.array(0)
                for k in range(np.hsplit(i.get_population().get_f(),2)[1].size):
                    _ = np.append(_, y[k] - i.get_population().get_x()[k][3]/1000)
    
                y = _[1:]
                del _           

        else: 
            x = np.array(0)
            for k in range(i.get_population().get_f().size):
                x = np.append(x, i.get_population().get_x()[k][i.get_population().get_x()[0].size - 1] / 365.25)    #years
            x = x[1:]
            y = i.get_population().get_f()/1000  #km/s delta v

            if add_vinf_arr==True:
                _ = np.array(0)
                for k in range((i.get_population().get_f().size)):
                    _ = np.append(_, y[k] - i.get_population().get_x()[k][3]/1000)
    
                y = _[1:]
                del _

        plt.figure(dpi,figsize=fig_size)
        plt.scatter(x, y, marker_size, c=y)
            
        if i.get_population().get_f()[0].size >  1:
            plt.title("Multi-Objective Optimized Solutions")
        else:
            plt.title("All Evolved Solutions of Island %i" %j)
        
        plt.xlabel("Time (years)")
        plt.ylabel("Delta V (km/s)")
        cbar= plt.colorbar()
        cbar.set_label("Delta V", labelpad=+3)
        plt.show()
        j = j+1   

# Keplerian Elements based ephemeris (epoch - 18/5/21)
earth_kep = pk.planet.keplerian((pk.epoch(2459352.5000,"jd")),(1.000245016541887 * pk.AU, 1.628852515836357E-02, 2.243694217384195E-03 * pk.DEG2RAD, 1.521366133881882E+02 * pk.DEG2RAD, 3.105719022836347E+02 * pk.DEG2RAD, 1.329425101759365E+02 * pk.DEG2RAD), pk.MU_SUN, 398600.435436*10**9, 6378137, 6569481.11, 'earth')
jupiter_kep = pk.planet.keplerian((pk.epoch(2459352.5000,"jd")),(5.20298651968408 * pk.AU, 4.839916594951266E-02, 1.303649764177414 * pk.DEG2RAD, 1.005193377709736E+02 * pk.DEG2RAD, 2.734208449481092E+02 * pk.DEG2RAD, 3.091406637921353E+02 * pk.DEG2RAD), pk.MU_SUN, 126686531.900*10**9, 69911000, 349555000, 'jupiter')
saturn_kep = pk.planet.keplerian((pk.epoch(2459352.5000,"jd")),(9.580106002660925 * pk.AU,5.201695072986685E-02, 2.482745756395404 * pk.DEG2RAD, 1.135751142783913E+02 * pk.DEG2RAD, 3.362338402237527E+02 * pk.DEG2RAD, 2.215098393207157E+02 * pk.DEG2RAD), pk.MU_SUN, 37931206.159*10**9, 58232000, 145580000, 'saturn')

seq = [earth_kep, jupiter_kep, saturn_kep]

earth_dep_1 = pk.epoch_from_string("2025-01-01 12:00:00")
earth_dep_2 = pk.epoch_from_string("2051-01-01 12:00:00")

udp = pk.trajopt.mga_1dsm(
        seq=seq,
        t0=[earth_dep_1, earth_dep_2],
        tof_encoding="alpha",
        tof=[8.0 * 365.25, 10 * 365.25],
        vinf=[7.0, 12.0],
        add_vinf_dep= True,
        add_vinf_arr= True,
        multi_objective=False,
        orbit_insertion=True,
        rp_target = saturn_kep.radius * 1.705,
        e_target = 0.99,
    )

prob = pg.problem(udp)
print(prob)

# Using Islands to create an Archipelago and different algos are used in different islands.
algo_archi = pg.algorithm(pg.mbh(pg.de(gen=1000,F=0.8,CR=0.9,variant=2)))  #Monotonic Basin Hopping with Differential Evolution

algo_isl_1 = pg.algorithm(pg.sade(gen=1000))         #Self Adaptive Differntial Algo
algo_isl_2 = pg.algorithm(pg.pso_gen(gen=1000,omega=0.715,eta1=1.7,eta2=1.7,variant=5,neighb_type=2))         #Particle Swarm Algo
algo_isl_3 = pg.algorithm(pg.sga(gen=1000,cr=0.5,m=0.091))          #Simpe Genetic Algo
algo_isl_4 = pg.algorithm(pg.de(gen=1000,F=0.8,CR=0.9,variant=2))

archi = pg.archipelago(algo = algo_archi, prob = prob, t=pg.ring(n=4,w=1))
isl_1 = pg.island(algo=algo_isl_1, prob=prob, size=100)
isl_2 = pg.island(algo=algo_isl_2, prob=prob, size=100)
isl_3 = pg.island(algo=algo_isl_3, prob=prob, size=100)
isl_4 = pg.island(algo=algo_isl_4, prob=prob, size=200)

# Calling relevant functions to run the program consisting of Single Objective isl and archi together.   
isl = [isl_1,isl_2,isl_3,isl_4]
size = run_isl(isl, [archi], udp)
scatter_all_sols(isl,marker_size=41, add_vinf_arr=True)
run_archi([archi], udp, size)

# Uncomment this section, and comment out the above code till seq, after running the code first time.
seq_EJ = [earth_kep, jupiter_kep]

# Get these values from Earth-Jupiter-Saturn output
earth_dep_EJ_1 = pk.epoch_from_string("2031-Feb-24 22:00:00.561631")
earth_dep_EJ_2 = pk.epoch_from_string("2031-Feb-24 22:45:55.561631")

udp_EJ = pk.trajopt.mga_1dsm(
        seq=seq_EJ,
        t0=[earth_dep_EJ_1, earth_dep_EJ_2],
        tof_encoding="alpha",
        tof=[2.50 * 365.25, 2.531 * 365.25],    # Get these values from Earth-Jupiter-Saturn output
        vinf=[8.4966, 8.5000],                  # Get these values from Earth-Jupiter-Saturn output
        add_vinf_dep= True,
        add_vinf_arr= True,
        multi_objective=False,
        orbit_insertion=True,
        rp_target = jupiter_kep.radius * 1.06,
        e_target = 0.99,
    )

prob_EJ = pg.problem(udp_EJ)
print(prob_EJ)

# Using Islands to create an Archipelago and different algos are used in different islands.
algo_archi_EJ = pg.algorithm(pg.mbh(pg.de(gen=1000,F=0.8,CR=0.9,variant=2)))  #Monotonic Basin Hopping with Differential Evolution

algo_isl_EJ_1 = pg.algorithm(pg.sade(gen=1000))         #Self Adaptive Differntial Algo
algo_isl_EJ_2 = pg.algorithm(pg.pso_gen(gen=1000,omega=0.715,eta1=1.7,eta2=1.7,variant=5,neighb_type=2))         #Particle Swarm Algo
algo_isl_EJ_3 = pg.algorithm(pg.sga(gen=1000,cr=0.5,m=0.091))          #Simpe Genetic Algo
algo_isl_EJ_4 = pg.algorithm(pg.de(gen=1000,F=0.8,CR=0.9,variant=2))

archi_EJ = pg.archipelago(algo = algo_archi_EJ, prob = prob_EJ, t=pg.ring(n=4,w=1))
isl_EJ_1 = pg.island(algo=algo_isl_EJ_1, prob=prob_EJ, size=100)
isl_EJ_2 = pg.island(algo=algo_isl_EJ_2, prob=prob_EJ, size=100)
isl_EJ_3 = pg.island(algo=algo_isl_EJ_3, prob=prob_EJ, size=100)
isl_EJ_4 = pg.island(algo=algo_isl_EJ_4, prob=prob_EJ, size=200)

# Calling relevant functions to run the program consisting of Single Objective isl and archi together.   
isl_EJ = [isl_EJ_1,isl_EJ_2,isl_EJ_3,isl_EJ_4]
size_EJ = run_isl(isl_EJ, [archi_EJ], udp_EJ)
scatter_all_sols(isl_EJ, marker_size=41, add_vinf_arr=True)
run_archi([archi_EJ],udp_EJ, size_EJ)

