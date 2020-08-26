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


CS.setup.closeloop(7,1)                                                         # Adding the equation of the pressure equivalence between the streams 1 and 7 


CS.setup.parameterset(1,"T",300)                                                # Temperature of the first stream = 300 K
CS.setup.parameterset(5,"T",299)                                                # Temperature of the 5th stream = 299 K
CS.setup.parameterset(3,"T",60)                                                 # Temperature of the 3th stream = 60 K
CS.setup.parameterset(4,"T",70)                                                 # Temperature of the 4th stream = 70 K


CS.setup.initvaluesset( ([299, 100, 80, 90, 297, 350, 300], [3, 3, 1, 1, 1, 3, 3], [1,1,1,1,1,1,1])) # Setting the initial values of cariables in form of ([T1,T2,...T7],[P1,P2,...P7],[m1,m2,...,m7])



"""
Optimzer:
"""
def Loptimize(variable):
    CS.setup.parameterset(3,"P", variable)                                      # Setup of pressure as current value of a solver variable (equation of pressure balance is added into CS.mainvar.funclist)
    CS.solution("TS")
    del CS.mainvar.funclist[-1]                                                 # Delete the last equation of the pressure balance from the system
    COP = HX2.capacity / C1.PowerF()                                            # Calculation of cycle efficiency
    return (1/COP)

pressure = 2                                                                    # Initial value of pressure

# Optimization of turbine outlet pressure within bounds from 1 to 30 bar
# The module scipy.optimize is already imported in the solver, so scipy.optimize.minimize is called from the solver. It is possible to import it directly in the runninng file.
solution = CS.opt.minimize(Loptimize, pressure, bounds=((1,30),),method='L-BFGS-B',options={'eps':0.0001}) 

print('Coefficient of performance:', (1/solution.fun) )
print('Pressure at the turbine outlet, bar:', solution.x)




### REMOVE IF PACKAGE IS INSTALLED
sys.path.remove(path)
### REMOVE IF PACKAGE IS INSTALLED