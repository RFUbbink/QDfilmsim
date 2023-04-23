# -*- coding: utf-8 -*-
"""
Created on Wed Apr 19 14:28:07 2023

@author: reinoutubbink
"""

#Rudimentary DOS drawer. I did not feel like writing it into a proper script with a config file etc. etc. 
#Just adapt parameters in the script itself if you want to use it or write your own DOS.
#The ZnO and CdSe parameters reproduce the DOS function used in the paper. 

import numpy as np
import math
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit

def Gauss(x,mu,sigma):
    return 1/math.sqrt(math.pi*2)/sigma*math.exp(-(x-mu)**2/2/sigma**2)

def threeOverTwo(x,a):
    return a*x**(3/2) 

def DOS_function(energy,energyLevels=[-4.35,-4.28,-3.7,-0.6],sigmas=[0,0.07,0.16,0.05],prefactors=[45,0.3,3,0]): 
    '''
    Return the DOS at a given energy.
    Can be tuned by adjusting the Energy levels, sigmas and prefactors of the sqrt contribution and 3 guassian contributions. 
    EnergyLevels: relative to vacuum level. 
    EnergyLevels and sigmas: in eV
    '''
    DOS = np.zeros(len(energy))
    for i,E in enumerate(energy):
        if E < energyLevels[0]:
            continue #To avoid imaginary numbers from the square root mostly. 
        else:
            energySum = 0
            energySum += prefactors[0]*math.sqrt(E-energyLevels[0]) #Add the sqrt contribution
            for j in range(1,4): #Add the guassian contributions
                energySum += prefactors[j]*Gauss(E,energyLevels[j],sigmas[j])
        DOS[i] = energySum
    return DOS

if __name__ == "__main__":
    DOS = [-4.6,0.01,0] 
    DOSint = [0.0] 
    N0 = 1.2e25 #Number of QDs per m-3 (was optimized based on experimental results)
    '''
    DOS[0] is always the starting energy relative to the vacuum level in eV!!
    DOS[1] is always the energy step size
    DOS[2] is the fit parameter for electron energy level correction (see below)
    DOSint is the integral over the DOS, which will be used later for fitting the fit parameter just mentioned
    '''
    energy = np.arange(DOS[0],DOS[0]+1.9,DOS[1]) #create the energy array. When scanning up to a voltage of -1 V, going up to 1.9 eV is more than enough to properly determine electron concentration
    
    '''
    For ZnO DOS:
    energyLevels=[-4.35,-4.28,-3.7,-2],sigmas=[0,0.07,0.16,0.1],prefactors=[45,0.3,3,0] 
    for CdSe DOS:
    energyLevels=[-10,-4.3,-4,-3.7],sigmas=[0,0.07,0.1,0.1],prefactors=[0,2,6,15] 
    '''
    DOS = np.append(DOS,N0*DOS_function(energy,energyLevels=[-4.35,-4.28,-3.7,-0.6],sigmas=[0,0.07,0.16,0.05],prefactors=[45,0.3,3,0]))
    for i in range(len(DOS)-3):
        DOSint.append(DOSint[i]+DOS[i+3]*DOS[1])
        
    '''
    Fitting the electron DOS correction factor
    This one is a doozy and it does improve the fit with experiment noticeably but not by that much
    It works as follows: 
    1) We fill up the DOS with more and more electrons. If there are more electrons in a certain part of the QD film, they will be in higher energy levels. 
    2) So the Fermi level in parts with more electrons lies higher, giving rise to additional force (besides normal diffusion) that pushes the electrons towards regions with less concentration.
    3) We want to know the difference in Fermi level between each cell, so that we can easily caculate the "total" energy level of electrons (=Fermi level + electrostic potential)
    4) We then use that total energy level instead of just the electrostatic potential level in a cell for the drift-diffusion equations
    5) Hope you are still fowllowing this. 
    6) So we want to know the Fermi level in each cell, which we can get by integrating the DOS, but this is computationally expensive. If we need to do this for every cell every timestep, we slow down the simulation A LOT
    7) So instead we hope to find a function that can be used to calculate the Fermi level directly from the electron concentration.
    8) Biggest assumption: DOS function is sqrt with energy. This is almost true for ZnO, not so much for CdSe maybe but it still fits rather OK.
    9) So then doing the math the integral of the DOS function should be a function of E^(3/2)
    10) So electron concentration scales with E^(3/2)
    11) And we can fit the reverse function to the electron concentration to obtain the Fermi energy
    12) Fermi level ~ electron concentration^(2/3)
    13) What we really want is the difference between 2 energy levels though, so in the simulator, we use the differential function dE/dn ~ n^(-1/3))
    14) So you will find in the simulator that dE = fit_factor*dn/n^(1/3). Which is really complicated but it is also the fastest calculation possible.
    15) If you feel this is too complicated, you can easily turn this functionality off by adjusting some things in the cpp source code. 
    '''
    
    energyForFitting = np.arange(0,1,DOS[1])
    startingIndexDOS = 3
    while DOS[startingIndexDOS] < 1:
        startingIndexDOS += 1
    popt,pcov = curve_fit(threeOverTwo,energyForFitting[:60],DOSint[startingIndexDOS:startingIndexDOS+60],[0.8e26]) #Fit the first 0.6 eV of the integrated DOS, using the energy relative to the LUMO
    #popt[0] = N0*45*2/3 #Alternatively, you can forego the fitting and use the real paramaters (this is wat is fitted really). This should give the same results if you have a pure sqrt DOS
    DOS[2] = (1/popt[0])**(2/3)*2/3 #Convert to proper unit for use in the simulator
    np.savetxt("DOS.csv", DOS, delimiter='\n')
    
    
    #Just plotting the DOS for quick check
    fig1,ax1 = plt.subplots(1,1)
    ax1.set_title("Here is your DOS")
    plt.plot(energy,DOS[3:],linewidth=2,color='k')
    plt.xlim([-4.6,-3.6])
    ax1.set_xlabel('Energy vs vacuum (eV)', size = '15')                              
    ax1.set_ylabel('DOS (m-3)', size = '15')
    ax1.tick_params(which='major',length=8, width=2,labelsize='15')
    plt.show()
    
    
    #And plotting the fit of the electron correction parameter
    fig2,ax2 = plt.subplots(1,1)
    ax2.set_title("Parameter fit for checking")
    plt.plot(energyForFitting[:60],DOSint[startingIndexDOS:startingIndexDOS+60],'k')
    plt.plot(energyForFitting[:60],threeOverTwo(energyForFitting[:60],*popt),'b')
    ax2.set_xlabel("Energy above LUMO", size = '15')                              
    ax2.set_ylabel("Electron concentration (m-3)", size = '15')
    plt.legend(["Data","Fit"])
    
    