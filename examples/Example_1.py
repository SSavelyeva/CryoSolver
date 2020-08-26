# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:25:46 2020

CryoSolver package - Cycle simulation example
@author: Sofiya Savelyeva, Steffen Klöppel
TU Dresden 2020
"""
"""
.........................................................................

Simple reverse Brayton cycle with 1 counter-current heat exchanger, 1 turbine, 1 load heat exchanger, 1 compressor and 1 aftercooler

"""
### REMOVE IF PACKAGE IS INSTALLED
import sys
path = '../'
if (path not in sys.path):
    sys.path.append(path)
### REMOVE IF PACKAGE IS INSTALLED

import cryosolver.solver as CS                                                  # Solver import


CS.start()                                                                      # Reset of solver global variables
CS.setup.fluidsetup("Helium.FLD")                                               # Setup of the working fluid 

HX1=CS.Heat_exchanger([CS.Cell(1,2,0.2),CS.Cell(4,5,0.1)])                      # Definition of the inner heat exchanger with pressure drops of 0.2 and 0.1 bar for high and low pressure sides
T1=CS.Turbine(2,3,0.85)                                                         # Definition of the turbine with isentropic efficiency of 85 %
HX2=CS.Simple_heat_exchanger(3,4,0.1,5)                                           # Definition of a load heat exchanger with a pressure drop of 0.1 bar and heat load of 5 kW
C1=CS.Compressor(5,6,0.8)                                                       # Definition of a compressor with an isentropic efficiency of 80 %
HX3=CS.Simple_heat_exchanger(6,7,0.1, Tout=300)                                   # Definition of a compressor aftercooler with a pressure drop of 0.1 bar and outlet temperature of 300 K


CS.setup.closeloop(1,7)                                                         # Adding the equation of the pressure equivalence between the streams 1 and 7 


CS.setup.parameterset(1,"T",300)                                                # Temperature of the first stream = 300 K
CS.setup.parameterset(5,"T",299)                                                # Temperature of the 5th stream = 299 K
CS.setup.parameterset(3,"P",2)                                                  # Pressure of the 3th stream = 2 bar
CS.setup.parameterset(3,"T",60)                                                 # Temperature of the 3th stream = 60 K
CS.setup.parameterset(4,"T",70)                                                 # Temperature of the 4th stream = 70 K

CS.setup.initvaluesset( ([299, 100, 80, 90, 297, 350, 300], [3, 3, 1, 1, 1, 3, 3], [1,1,1,1,1,1,1])) # Setting the initial values of cariables in form of ([T1,T2,...T7],[P1,P2,...P7],[m1,m2,...,m7])

s = CS.solution("TS")                                                              # Solving the equation by means of scipy.optimize.root  

CS.output.txtconsole(s)                                                         # Console output in the text form

CS.output.exergylosses()                                                        # Show distribution of exergy losses


"""
User output examples:
"""
print('Temperature of the 2st stream in K:', CS.mainvar.streamlist[1].T)        # Here index of array is used, which is i-1 for the i-stream
print('Compressor power, kW: ', C1.power)
print("Heat exchanger, capacity of high pressure stream, kW:", HX1.cells[0].capacity)  # Here index 0 states for the capacity of the high pressure side
print("Heat exchanger, capacity of low pressure stream, kW:", HX1.cells[1].capacity)   # Here index 1 states for the capacity of the low pressure side
print('Aftercooler capacity, kW: ', HX3.capacity)

print(HX1.HXdef(Para='mindt', savedata=False, plot=True))                       # Printing the heat exchanger Q-T diagram; savedata=True saves the .txt file with heat load and temperature values in the examples folder

CS.output.excel()                                                               # Saving a new excel protocol file in the examples folder

### REMOVE IF PACKAGE IS INSTALLED
sys.path.remove(path)
### REMOVE IF PACKAGE IS INSTALLED
