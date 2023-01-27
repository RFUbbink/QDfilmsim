# -*- coding: utf-8 -*-
"""
Created on Wed Dec 23 12:20:09 2020

@author: reinoutubbink
"""

import numpy as np
import matplotlib.pyplot as plt

dt = 1e-10      #time step size

Vstart = -0.1      #bias applied to the working electrode vs. reference
Vstop = -0.6
Sweep = 10000   #sweep rate in Vs-1
incr = 0.01     #sweep increment V
tstep = int(incr/Sweep/dt)#time after which the voltage is increased
num = int(abs(Vstop-Vstart)/incr) #number of increments that need to be taken (one way)
S = 2           #SCAN TYPE [1]linear sweep [2]cyclic voltammetry

#%% read C++ output file and analyse

JStore = np.loadtxt('JStore.csv', delimiter=',')

plt.figure('ncurrent')
plt.plot(JStore)

#plt.plot(np.arange(Vstart, Vstop, -incr),JStore[:num],'b')
#plt.plot(np.arange(Vstop+incr, Vstart+incr, incr),JStore[num:],'b')
#plt.xlabel('WE bias (V)', size = '15')
#plt.ylabel('Current density (A/cm^2)', size = '15')
#plt.tick_params(labelsize='15')
#plt.show