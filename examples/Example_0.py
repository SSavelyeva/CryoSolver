# -*- coding: utf-8 -*-
"""
CryoSolver package - Simulation example
@author: Sofiya Savelyeva, Steffen Klöppel
TU Dresden 2020
"""
"""
.........................................................................

Calculation of a turbine

"""
### REMOVE IF PACKAGE IS INSTALLED
import sys
path = '../'
if (path not in sys.path):
    sys.path.append(path)
### REMOVE IF PACKAGE IS INSTALLED


import cryosolver.solver as CS                                                  # Import of the solver
CS.start()                                                                      # Creation of global variables

CS.setup.fluidsetup("Nitrogen.FLD","Helium.FLD")                                # Fluid definition

CS.Turbine(1,2,0.85,power=700)                                                  # Definition of the turbine with 85 % isentropic efficiency and power of 700 kW


CS.setup.parameterset(1,"Z",[0.6,0.4])                                          # Setup of the composition of the inlet stream
CS.setup.parameterset(2,"T",100)                                                # Fixing the outlet temperature of the turbine as 100 K
CS.setup.parameterset(2,"P",1)                                                  # Fixing the outlet pressure of the turbine as 1 bar


def userfunction():                                                             # User function to fix the temperature drop through the turbine as 100 
    return CS.mainvar.streamlist[0].T - CS.mainvar.streamlist[1].T - 100   
CS.mainvar.funclist.append(userfunction)


CS.setup.initvaluesset( ([300,200], [5,1], [1,1], [[0.78],[0.78]])  )           # Initial values for the solver as array: ([T1,T2],[P1,P2],[m1,m2],[[z_nitrogen],[z_nitrogen]])

s = CS.solution("TP")                                                           # Finding the solution based on the temperature-pressure fluid properties calculation

CS.output.txtconsole(s)                                                         # Output of the solution in the text form on the console




### REMOVE IF PACKAGE IS INSTALLED
sys.path.remove(path)
### REMOVE IF PACKAGE IS INSTALLED