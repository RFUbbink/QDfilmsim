

#!/usr/bin/env pypy
# -*- coding: utf-8 -*-
"""
Created on Mon Oct 19 16:06:26 2020

@author: Gebruiker
"""

import numpy as np
import math
import matplotlib.pyplot as plt

#Define parameters
#spacetime
L = 2100e-9      #distance between WE and CE 
Ref = 1400e-9    #reference electrde position
Fth = 700e-9    #QD film thickness 
N = 90          #amount of cells
H = 1+math.ceil(Fth/L*N)#film/solution interface
R = int(Ref/L*N)#position of reference electrode
dx = L/(N-1)    #cell width (diameter of QD)
Awe = 1e-4      #area of the working electrode m2
Ace = 1e-4      #area of the counterelectrode m2
refFactor = Ref/L #for faster conversion in the Potential function
dt = 5e-9    #time step size

#scan charactierstics
Vstart = -0.1   #bias applied to the working electrode vs. reference
Vstop = -0.8
Sweep = 1000   #sweep rate in Vs-1
incr = 0.01     #sweep increment V
tstep = int(incr/Sweep/dt)#time after which the voltage is increased
num = abs(Vstop-Vstart)/incr #number of increments that need to be taken (one way)
S = 2           #SCAN TYPE [1]linear sweep [2]cyclic voltammetry


#mobilities
mobn = 1e-9    #electron mobility m2V-1s-1
mobaf = 0.5e-11    #anion mobility in the QD film m2V-1s-1
mobcf = 3e-12     #cation mobility in the QD film m2V-1s-1
mobas = 1e-9   #anion mobility in solution m2V-1s-1
mobcs = 1e-9    #cation mobility in solution m2V-1s-1

#misc
N0 = 5e25       #number of QD's m-3
c0 = 6e25    #salt concentration m-3
T = 300         #Temperature
epsrf = 10       #relative dielectric constant in quantum dot film
epsrs = 37      #relative dielectric constant of solution

#band characteristics
LUMO = -4.3         #eV relative to vacuum
nELWF = -4.7        #eV relative to vacuum
dE = 0.01           #eV for integration resolution
Emax = LUMO + 3     #eV upper bound for DOS intgration
intSteps = int((Emax-LUMO)/dE)
phin = LUMO - nELWF  #injection barriers = difference between band edges and work function of respective electrode material
Vce = Vstart*Awe/Ace*epsrf/epsrs #bias applied to counterelectrode vs. reference
Vb = Vce + Vstart        #total applied bias between the working electrode and reference

#natural constants
k = 1.3806e-23      #m2kgs-2K-1
eps0 = 8.8542e-12   #m-3kg-1s-4A2
h = 6.6261e-34      #m2kgs-1
q = 1.6022e-19      #C
m0 = 9.1094e-31     #kg

#complex constants constants for speed
eConvert = k*T/q
eConvertx = eConvert/dx
curConn = mobn*dt/dx
curConcf = mobcf*dt/dx
curConaf = mobaf*dt/dx
curConcs = mobcs*dt/dx
curConas= mobas*dt/dx
PotAf = -q/eps0/epsrf
PotAs = -q/eps0/epsrs


#creating initial electrochemical CELL
CELL = np.zeros((3,N)) 
#[0]electrons [1]cations [2]anions
potCELL = np.zeros((3,N))
#[0] ddV [1]electric field [2]electrostatic potential 
curCELL = np.zeros((3,(N+1)))
#[0]electrons [1]cations [2]anions
CELL[1][1:H],CELL[2][1:H] = c0*0.2, c0*0.2                         #set the starting concentrations 
CELL[1][H:(N-1)],CELL[2][H:(N-1)] = c0, c0                      #concentration in film is half that in solution because QDs take up space
potCELL[2][0] = Vstart                                          #setting a boundary condition for the ES potential                                  

JStore = []                                                     #to store the current

Jstore2 = []
#defining functions
def Curn(conL,conR,curCon,EF): #outputs amount of negative charges moving to right cell per time step per m3

    return (-conL*EF + eConvertx*(conL-conR))*curCon     #based on drift/diff equation adapted to proper units for this code


def Curp(conL,conR,curCon,EF): #outputs amount of positive charges moving to left cell per time step per m3
    return (conR*EF + eConvertx*(conL-conR))*curCon        #based on drift/diff equation adapted to proper units for this code


def injectBM():   #Modified Boltzman model for injcetion of electrons 
    return N0*math.exp((potCELL[2][1]-potCELL[2][0]-phin)/eConvert) 


def calcPot(Vb):                # Calculate Pot based on concentrations
    for x in range(1,H):
        potCELL[0][x] = PotAf*(-CELL[0][x]+CELL[1][x]-CELL[2][x])  #Poisson equation
    for x in range(H,(N-1)):
        potCELL[0][x] = PotAs*(-CELL[0][x]+CELL[1][x]-CELL[2][x])
    Vref = 1
    while abs(Vref) > 0.0001:           #this boundary condition makes sure that the potential at the reference is 0
        E0 = Vb/L                       #try initial EF boundary condition. Assume 'the rest' of the potential (i.e. everything not due to accumulated charges) drops linearly over the active layer. 
        diffE = 1
        while abs(diffE) > 0.00000001:       
            potCELL[1][0] = E0                #Build Electric field based on second space derivatives
            for i in range(N-2):
                potCELL[1][i+1] = potCELL[1][i] - potCELL[0][i+1]*dx
            Eint = sum(potCELL[1])*dx
            diffE = (Vb-Eint)/Vb       #implementation of boundary condition: Total potential drop = Vbias = integral of EF to dx
            E0 = E0*(1+diffE)
    
        for i in range(N-1):                 #Build potential based on first space derivative (electric field)
            potCELL[2][i+1] = potCELL[2][i] - potCELL[1][i]*dx
        Vref = potCELL[2][R]
        Vb += (Vref/refFactor) 
#        print('Vref = ',Vref)
    return Vb

def calcCurrents():
    for i in range(1,H): #charge transport electrons
        curCELL[0][i] = Curn(CELL[0][i-1],CELL[0][i],curConn,potCELL[1][i-1])            #calculate electron currents between each cell
    for i in range(2,H):
        curCELL[1][i] = Curp(CELL[1][i-1],CELL[1][i],curConcf,potCELL[1][i-1])           #calculate anion currents between each cell
        curCELL[2][i] = Curn(CELL[2][i-1],CELL[2][i],curConaf,potCELL[1][i-1])           #similar for cations
    for i in range((H+1),(N-1)): #charge transport ions solution
        curCELL[1][i] = Curp(CELL[1][i-1],CELL[1][i],curConcs,potCELL[1][i-1])            
        curCELL[2][i] = Curn(CELL[2][i-1],CELL[2][i],curConas,potCELL[1][i-1])
#    curCELL[1][H] = (-math.sqrt(CELL[1][H-1]*0.2*CELL[1][H])*potCELL[1][H-1] + eConvertx*(CELL[1][H-1]-CELL[1][H]))*curConcf
#    curCELL[2][H] = (math.sqrt(CELL[2][H-1]*0.2*CELL[2][H])*potCELL[1][H-1] + eConvertx*(CELL[2][H-1]-CELL[2][H]))*curConaf        
    curCELL[1][H] = Curp(CELL[1][H-1],CELL[1][H]*0.2,curConcf,potCELL[1][H-1]) 
    curCELL[2][H] = Curn(CELL[2][H-1],CELL[2][H]*0.2,curConaf,potCELL[1][H-1]) 
#    curCELL[1][H] = -math.exp(-(0.15 - potCELL[2][H] + potCELL[2][H-1])*q/k/T)*CELL[1][H] + CELL[1][H-1];
#    curCELL[2][H] = -math.exp(-(0.15 + potCELL[2][H] - potCELL[2][H-1])*q/k/T)*CELL[2][H] + CELL[2][H-1];
   
def updateConc():
    for h in range(1,H):
        CELL[0][h] += curCELL[0][h] -curCELL[0][h+1]                                        #update concentrations using currents
    for h in range(1,(N-1)):
        CELL[1][h] += curCELL[1][h] -curCELL[1][h+1]                                        #update ion concentrations 
        CELL[2][h] += curCELL[2][h] -curCELL[2][h+1]
    CELL[1][R],CELL[2][R] = c0,c0


#%% Run code
#CELL = np.loadtxt('Outputpy.csv', delimiter=',')               #If you want to load an earlier experiment
steps = int(num*S)
#Vb = -0.402095747629599
for V in range(steps):
    for t in range(tstep):
        Vb = calcPot(Vb)                                                    #calculate EP and EF. Also updates Vb to speed up next timestep         
        CELL[0][0] = injectBM()                                           #injection of electrons 
        calcCurrents()
        updateConc()
    JStore.append(-curCELL[0][2]*q*dx/dt/10000)                     #after each increment, record the current density (converted to A/cm2)                                 
    print(V)
    if V < steps/2:                                                         #inreases/decreases the applied bias with the increment after a certain amount of time
        Vb -= incr
        potCELL[2][0] -= incr
    else:
        Vb += incr
        potCELL[2][0] += incr
#np.savetxt('Outputpy.csv', JStore, delimiter=',')
        
#%%
#print(CELL)  
dis = np.linspace(0,(L*1e9),N)
plt.close('all')
plt.figure('Carrier concentrations')               #Plots ion and electron/holes concentrations throughout the device
plt.plot(dis,CELL[0][0:N],'b')
plt.plot(dis[1:(N-10)],CELL[1][1:(N-10)],'r--')
plt.plot(dis[1:(N-10)],CELL[2][1:(N-10)],'b--')
plt.xlabel('Distance (nm)', size = '15')
plt.ylabel('Carrier concentrations (m-3)', size = '15')
plt.legend(['electrons','cations','anions'])
plt.figure('Electrostatic potential')
plt.plot(dis,potCELL[2])
plt.xlabel('Distance (nm)', size = '15')
plt.ylabel('Electrostatic potential (V/m)', size = '15')
plt.tick_params(labelsize='15')
plt.figure('ncurrent')
plt.plot(np.arange(Vstart, Vstop, -incr),JStore[:int(num)],'b')
plt.plot(np.arange(Vstop+incr, Vstart+incr, incr),JStore[int(num):],'b')
plt.xlabel('WE bias (V)', size = '15')
plt.ylabel('Current density (A/cm^2)', size = '15')
plt.tick_params(labelsize='15')
plt.show

#np.savetxt('QDfilmleakQS.csv', CELL, delimiter=',')
