# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:25:46 2020

CryoSolver package - Cycle simulation example
@author: Sofiya Savelyeva, Steffen Klöppel
TU Dresden 2020
"""
"""
.........................................................................

Claude cycle

"""
### REMOVE IF PACKAGE IS INSTALLED
import sys
path = '../'
if (path not in sys.path):
    sys.path.append(path)
### REMOVE IF PACKAGE IS INSTALLED


import cryosolver.solver as CS                                                  # Solver import

CS.start()                                                                      # Reset of solver global variables
CS.setup.fluidsetup("Nitrogen.FLD","Oxygen.FLD")                                # Setup of the working fluid

HX1=CS.Heat_exchanger([CS.Cell(1,2,0.1),CS.Cell(10,11,0.1)])                    # Definition of the inner heat exchanger with 2 flows and 0.1 bar pressure drop of each side
HX2=CS.Heat_exchanger([CS.Cell(3,4,0.1),CS.Cell(9,10,0.1)])                     # Definition of the inner heat exchanger with 2 flows and 0.1 bar pressure drop of each side
HX3=CS.Heat_exchanger([CS.Cell(4,5,0.1),CS.Cell(7,8,0.1)])                      # Definition of the inner heat exchanger with 2 flows and 0.1 bar pressure drop of each side
Sep=CS.Separator(2,15,3)                                                        # Definition of a flow separator with inlet with index 2 and outlets 15 and 3
Mix1=CS.Mixer(9,8,16)                                                           # Definition of a mixer with an outlet with index 9 and inlets 8 and 16
T1=CS.Turbine(15,16,eta_is=0.85)                                                # Definition of a turbine with isentropic efficiency of 85 %
V=CS.Valve(5,6)                                                                 # Definition of a throttle valve
C1=CS.Compressor(12,13,eta_is=0.9)                                              # Definition of a compressor with 90 % isentropic efficiency
AC1=CS.Simple_heat_exchanger(13,14,0,Tout=300)                                    # Definition of a simple after-cooler with outlet temperature of 300 K
PS=CS.PhaseSeparator(6,7,17)                                                    # Definition of a phase separator with inlet 6, vapor outlet 7 and liquid outlet 17
Mix2=CS.Mixer(12,11,18)                                                         # Definition of a mixer with an outlet with index 12 and inlets 11 and 18

CS.setup.closeloop(1,13)                                                        # Adding the equation of the pressure equivalence between the streams 1 and 13 

CS.setup.parameterset(1, "Z",[0.78,0.22])                                       # Setting the volumetric mixture composition (for a mixture with N fluids only (N-1) mole fractions will me used, so it is possible to write [0.78] here instead of [0.78,0.22])
CS.setup.parameterset(1, "T", 300)                                              # Setting the temperature of the first stream as 300 K
CS.setup.parameterset(1, "M", 1)                                                # Setting the mass flow of the first stream as 1 kg/s
CS.setup.parameterset(15, "M", 0.9)                                             # Setting the mass flow through the turbine as 0.9 kg/s
CS.setup.parameterset(1, "P", 60)                                               # Setting the pressure of the first stream as 60 bar
CS.setup.parameterset(6, "P", 10)                                               # Setting the pressure at the turbine outlet as 10 bar
CS.setup.parameterset(2, "T", 220)                                              # Setting the temperature of the 2nd stream as 220 K
CS.setup.parameterset(5, "T", 140)                                              # Setting the temperature of the 5nd stream as 140 K
CS.setup.parameterset(12, "Z",[0.78,0.22])                                      # Setting the volumetric mixture composition at the inlet of compressor after the mixer

def func():                                                                     # Adding a custom function to define equality between temperatures of stream 8 and 16 (here indernal indexes of solver are used in functions, which need to be indicated as 8 - 1, 16 - 1)
    return CS.mainvar.streamlist[7].T - CS.mainvar.streamlist[15].T
CS.mainvar.funclist.append(func)

def func2():                                                                    # Adding a custom function to define equality between temperatures of stream 18 and 11 (here indernal indexes of solver are used in functions, which need to be indicated as 18 - 1, 11 - 1)
    return CS.mainvar.streamlist[17].T - CS.mainvar.streamlist[10].T
CS.mainvar.funclist.append(func2)

def func3():                                                                    # Adding a custom function to define equality between mass flows of stream 18 and 17 (here indernal indexes of solver are used in functions, which need to be indicated as 18 - 1, 17 - 1)
    return CS.mainvar.streamlist[17].m - CS.mainvar.streamlist[16].m
CS.mainvar.funclist.append(func3)




CS.setup.initvaluesset( ([300,220,220,160,140,120,120,130,130,138,299,299,400,300,220,138,120,299], [60,60,60,60,60,10,10,10,10,10,10,10,60,60,60,10,10,10],[1,1,0.1,0.1,0.1,0.1,0.06,0.06,0.96,0.96,0.96,0.96,1,1,0.9,0.9,0.04,0.06], [[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7],[0.7]]) ) # Setting the initial values of cariables in form of ([T1,T2,...T7],[P1,P2,...P7],[m1,m2,...,m7], [[Z1(N2)],[Z1(N2)],...[Z7(N2)]])

s=CS.solution("TS")                                                             # Solving the equation by means of scipy.optimize.root

CS.output.txtconsole(s)                                                         # Console output in the text form

CS.output.excel()                                                               # Saving an excel protocol

### REMOVE IF PACKAGE IS INSTALLED
sys.path.remove(path)
### REMOVE IF PACKAGE IS INSTALLED
