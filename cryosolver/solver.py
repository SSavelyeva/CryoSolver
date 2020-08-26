# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:25:46 2020

CryoSolver package
@author: Sofiya Savelyeva, Steffen Klöppel
TU Dresden 2020
"""

"""
.........................................................................
solver.py
The library of functions and classes to setup and simulate the cycle

"""
from sys import exit
import scipy.optimize as opt
import numpy as np
from math import exp, log
import cryosolver.refpropfunc as RF

"""
.........................................................................
Simulation setup and main variable definition
"""

class setup:
    global mainvar
        
    def fluidsetup(*fluids):
        """
        Function to setup the fluid throught the refpropfunc.py from a given number of pure fluids.
        *fluids -- multiple arguments of pure fluid names ("name.FLD","name2.FLD",...)
        """
        fluidlist=[]
        for fluid in fluids:
            fluidlist.append(fluid)                                             # Definition of a list for a given number of fluid names
        mainvar.fluid = RF.fluid(fluidlist)                                     # Fluid setup by means of class fluid in refpropfunc.py
        return          
    
    class mainvarsetup:  
        """
        Initial definition of the global "mainvar" class variable required for simulations.
        """
        def __init__(self,exergy0,exergycalculation,fixcomposition):
            self.funclist = []                                                  # List of bound methods of functions to collect the system of equations
            self.streamlist = []                                                # List to define all the streams of the cycle and their properties
            self.streamnumber = None                                            # Variable to store the number of defined streams (is equal to len(streamlist))            
            self.fluid = None                                                   # Variable to store the information on the used fluid (has a type of RF.fluid class)            
            self.guess = None                                                   # Array of initial values of temperatures, pressures, mass flows and compositions for the solver            
            self.componentslist = []                                            # List of class objects to store information about defined cycle components
            self.streaminletlist = []                                           # Utility list collecting all the streams entering the cycle components - for self-check
            self.streamoutletlist = []                                          # Utility list collecting all the outgoing streams of the cycle components - for self-check
            self.stype = None                                                   # String variable to setup the calculation type: "TS" for temperature-entropy based and "TP" for temperature-pressure based calcilations
            self.exergyzerostate = exergy0                                      # Zero state to calculate exergy (array in form [Temperture, pressure], by default [300,1])
            self.exergyflag = exergycalculation                                 # Boolean variable to set whether the exergy and exergy losses should be calculated (True by defualt)
            self.setupflag = False                                              # Boolean variable to set whether the setup.solversetup has been made
            self.fixcomposition = fixcomposition                                # Optional: only if the composition is known for all the cycle streams it can be specified in advance to reduce the number of system equations as array of compositions [[Z11,Z12,...],[Z21,Z22,...],...]
            return
        
    class parameterset: 
        """
        Class to setup known variables. It adds equations in form of (variable - value = 0) into the mainvar.funclist.
            stream_number - index of a cycle stream
            parameter - string of variable to be fixed ("T" for temperature, "P" for pressure, "M" ofr mass flow or "Z" for composition)
            value - corresponding value (float or array)
        """
        def __init__(self,stream_number, parameter, value):
            self.stream_number = stream_number - 1                              # Within the solver all indexes are staring with 0 instead of 1
            self.parameter = parameter
            self.value = value
            if self.parameter is "T":
                mainvar.funclist.append(self.functionT)                         # Add function (T-value) into mainvar.funclist
            elif self.parameter is "P":
                mainvar.funclist.append(self.functionp)                         # Add function (P-value) into mainvar.funclist
            elif self.parameter is "M":
                mainvar.funclist.append(self.functionm)                         # Add function (M-value) into mainvar.funclist
            elif self.parameter is "Z":
                """
                Within the solver the mole fractions are calculated for all the mixture components except the last one. 
                The mole fraction of a last mixture component is defined as (1-sum(Z[0:(mainvar.fluid.componentsnumber-1)])). 
                The composition of a pure fluid ([1.]) is set authomatically and is excluded from the system of equations.
                """
                if mainvar.fixcomposition is not None:
                    print("Composition of each stream will be specified through mainvar.fixcomposition")
                else:
                    for i in range(mainvar.fluid.componentsnumber - 1):                
                        setup.parameterset.funcinitialize(i, self.stream_number, self.value[i])      # Add function (Z[i]-value) into mainvar.funclist, where i is an index of mixture component
            else:
                print("Unknown parameter specified")
                exit()
            return
        def functionT(self):  
            return mainvar.streamlist[self.stream_number].T - self.value
        def functionp(self):        
            return mainvar.streamlist[self.stream_number].p - self.value
        def functionm(self):
            return mainvar.streamlist[self.stream_number].m - self.value
        class funcinitialize:
            def __init__(self, i, stream, value):
                self.i = i
                self.stream_number = stream
                self.value = value
                mainvar.funclist.append(self.functionz)
            def functionz(self):         
                return mainvar.streamlist[self.stream_number].z[self.i] - self.value

    class closeloop:
        """
        Function to setup equality of two streams to close the cycle loop.
        Adds the equation of pressure equality into mainvar.funclist. If necessary, the temperature equality can be included as well.
        Please note: if given streams are in 2-phase region, entropy equality should be setup instead of pressure equality.
        """
        def __init__(self,first_stream, second_stream):
            self.first_stream = first_stream - 1                                # Within the solver all indexes are staring with 0 instead of 1
            self.second_stream = second_stream - 1
            mainvar.funclist.append(self.MomentumB)
        def MomentumB(self):
            return mainvar.streamlist[self.first_stream].p - mainvar.streamlist[self.second_stream].p
        def TemperatureB(self):
            return mainvar.streamlist[self.first_stream].T - mainvar.streamlist[self.second_stream].T      # Additional function that can be added into the closeloop command: mainvar.funclist.append(self.TemperatureB)
    
    def initvaluesset(arrayofvalues = None):
        """
        Initialization of starting values for the solver
        Arrayofvalues should be given in form [[T1,T2,...],[P1,P2,...],[M1,M2,...],*[[Z11,Z12,..],[Z21,Z22],...]], 
        where T - temperature, P - pressure, M - mass flow, Z - mole fraction of mixture components except the last component (not required: for a pure fluid, if mainvar.fixedcomposition is True)
        """
        if arrayofvalues is None:
            exit("Please specify initial values for variables")
        else:
            mainvar.guess = arrayofvalues
        return 
    
    def packguess(guessarray):
        """
        Packing the mainvar.guess into a simple array of values depending on the solution type ("TS" or "TP") that will be used by a solver as starting values.
            guessarray has a form given in setup.initialvaluesset;
            guess has a form of simple 0-axis array: [T1,T2,...,P1 or S1,P2 or S2,...,M1,M2,...,*Z11,*Z12,..,*Z21,*Z22].
        """
        Tarray = guessarray[0]
        Parray = guessarray[1]
        Marray = guessarray[2]
        if mainvar.stype is "TS":
            if mainvar.fixcomposition is not None:                                  # If the composition is fixed by user, it will be used to calculate entropy
                Z = mainvar.fixcomposition
            elif mainvar.fluid.componentsnumber > 1:
                Zarray = guessarray[3]   
                Z = []                                                          # Z is the full composition including the mole fraction of the last mixture component (calculated here to define entropy)
                for i in range(mainvar.streamnumber):
                    Z.append(np.concatenate((Zarray[i], [1-sum(Zarray[i])])))
            else:
                Z = [[1.]*1]*mainvar.streamnumber 
            Sarray = np.array([None]*mainvar.streamnumber, dtype=float)         # Calculation of initial values for entropies for given initial temperatures, pressures and composition in case of temperature-entropy based simulations
            for i in range(mainvar.streamnumber):
                Sarray[i] = RF.S_TP(Z[i], Tarray[i], Parray[i])                 
            if (mainvar.fluid.componentsnumber==1) or (mainvar.fixcomposition is not None):
                guess = np.concatenate((Tarray, Sarray, Marray))
            else:
                guess = np.concatenate((Tarray, Sarray, Marray, np.concatenate((Zarray))))     
        else:
            if (mainvar.fluid.componentsnumber==1) or (mainvar.fixcomposition is not None):
                guess = np.concatenate((Tarray, Parray, Marray))
            else:
                Zarray = guessarray[3] 
                guess = np.concatenate((Tarray, Parray, Marray, np.concatenate((Zarray))))
        return guess
    
    
    def solversetup(): 
        """
        Initial setup of solver before simulation.
            1. It checks the number of assigned cycle streams and number of equations.
            Returns error if the total number of variables is not equal to the number of equations.
            2. In case of temperature-entropy based simulations it recalculates pressures into entropies in the mainvar.guess array.
        """
        def listelements(firstlist:list,secondlist:list):                       # Function to return elements of the first list that are not in the second list
            final_list = [] 
            for num in firstlist: 
                if num not in secondlist: 
                    final_list.append(num) 
            return final_list
        def removeduplicates(listvar:list):                                     # Function to remove the dupcilate numbers in two lists
            final_list = [] 
            for num in listvar: 
                if num not in final_list: 
                    final_list.append(num) 
            return final_list 
        joinstreams = sorted(removeduplicates((mainvar.streaminletlist + mainvar.streamoutletlist)))        # Total list of unique streams specified in the cicle 
        mainvar.streamnumber = len(joinstreams)                                 # Total number of specified streams stored in mainvar.streamlist
        print("Number of specified streams: ", mainvar.streamnumber)            # Program output for self-check. Please note: the open streams are shown independently on setup.closeloop() function initialisation.
        print("""Specified open streams:                                        
            inlets:""",  listelements(mainvar.streaminletlist, mainvar.streamoutletlist),                    
            """outlets:""", listelements(mainvar.streamoutletlist, mainvar.streaminletlist))
        if mainvar.fixcomposition is not None:                                      # Estimation of number of required variables (T,P/S,M,[Z1,Z2,...] for each stream): if the composition is fixed it is not used as variable in the solver
            required_variables = mainvar.streamnumber*3
        else:
            required_variables = mainvar.streamnumber*(3+(mainvar.fluid.componentsnumber-1))
        Unknownfunc = required_variables - len(mainvar.funclist)                # Difference between required variables and number of specified equations.
        if (Unknownfunc > 0) :
            print(Unknownfunc, " more function(s) need(s) to be specified.")
            exit("Programm stopped")
        elif (Unknownfunc < 0) :
            print("Too many functions are specified. Number of functions to delete:", np.abs(Unknownfunc) ) 
            exit("Programm stopped")
        if mainvar.guess is None:
            print("Please input initial values of variables!")
            exit("Programm stopped")
        mainvar.guess = setup.packguess(mainvar.guess)                          # Packing the mainvar.guess into a simple array of values depending on the solution type ("TS" or "TP")
        mainvar.setupflag = True
        return


def start(exergy0=[300,1], exergycalculation=True, fixcomposition = None):
    """
    Setup of initial parameters for the cycle simulation:
        1. Default definition of a "global mainvar"; 
            - "mainvar" is the global class object of variables required for the simulation. 
            - This is the main object for user to access stream properties within simulation.
            - It is initially defined through setup.mainvarsetup.__init__() (see setup.mainvarsetup for the list of available variables).      
        2. Setup of entropy calculation parameters  
            - exergy0 - zero state for exergy calculation in form of an array [Temperature,Pressure]. Default zero state is 300 K and 1 bar abs.
            - exergycalculation - boolean variable to setup weather the exergy losses are calculated (by default True).
        3. Optional: fixing the composition of all streams 
            - only if the composition is known for all the cycle streams it can be specified in advance to reduce the number of system equations, this will increase the calculation speed and accuracy
            - fixcomposition should have a for of array of compositions for each stream: [[Z11,Z12,...],[Z21,Z22,...],...]
            - shouldn't be used if a phase separator is specified
    """
    global mainvar
    mainvar = setup.mainvarsetup(exergy0, exergycalculation, fixcomposition)
    return
  

class stream:
    """
    The class to uniquely define any stream with its properties
        index - cycle stream index used within the solver (starting with 0) 
        T - temperature
        m - total mass flow
        z - full array of mole fractions for all the mixture components (or [1.] for a pure fluid)
        p - pressure
        s - entropy
        h - enthalpy
        ex - exergy related to the zero state defined in mainvar.exergyzerostate
    """
    global mainvar
    def __init__(self, T, m, z, index, p=None, s=None): 
        self.index = index 
        self.T = T            
        self.m = m                    
        self.z = z  
        if p is not None:
                self.p = p
                self.h = RF.H_TP(self.z, self.T, self.p)
                self.s = RF.S_TP(self.z, self.T, self.p)
        elif s is not None:
                self.s = s
                self.h = RF.H_TS(self.z, self.T, self.s)
                self.p = RF.P_TS(self.z, self.T, self.s)
        else:
            print("Error by stream calculation, please input pressure or entropy")
            exit("Programm stopped")
        if mainvar.exergyflag is False:
            self.ex = None
        else:
            self.ex = (self.h - RF.H_TP(self.z, mainvar.exergyzerostate[0], mainvar.exergyzerostate[1])) - mainvar.exergyzerostate[0]*(self.s-RF.S_TP(self.z, mainvar.exergyzerostate[0], mainvar.exergyzerostate[1]))
        return
   
    
    
"""
.........................................................................
Classes of components
"""    
class Compressor:
    """
    Input paramters:
        inlet_stream_number - index of inlet cycle stream starting with 1;
        outlet_stream_number - index of outlet cycle stream starting with 1;
        eta_is - isentropic efficiency (1 by default);
        pressure_ratio - p_outlet / p_inlet (optional);
        power - input power in kW (optional).
    Added equations:
        EnergyB - energy balance (h_out - h_in - (1/eta_is)*(h_outs - h_in))
        MassB - mass balance (m_in - m_out)        
        Optionally:
            CompositionB - composition equilibrium for each of mixture components except the last one (z_in[i] - z_out[i]), where i - component of the mixture; not defined for a pure fluid or if composition is fixed
            MomentumB - pressure balance if pressure ratio is given (p_out/p_in - pressure_ratio)
            PowerB - full energy balance if power is given (h_out*m_out - h_in*m_in - power)       
    PowerF, ExergyF, MomentumF - available functions to calculate actual power, exergy loss and pressure ratio
    """
    global mainvar 
    def __init__(self, inlet_stream_number, outlet_stream_number, eta_is = 1, pressure_ratio = None, power = None):
        self.name = "Compressor"
        self.eta_is = eta_is
        self.inlet = inlet_stream_number - 1                                    # Internal index of inlet stream starting with 0
        self.outlet = outlet_stream_number - 1                                  # Internal index of outlet stream starting with 0
        self.pressure_ratio = pressure_ratio
        self.power = power     
        self.exergy_loss = None
        mainvar.funclist.append(self.EnergyB)
        mainvar.funclist.append(self.MassB)
        if mainvar.fixcomposition is None:
            for i in range(mainvar.fluid.componentsnumber-1):
                Compressor.funcinitialize(i, self.inlet, self.outlet)           # Adding the composition equilibrium functions for each of mixture components except the last one
        if self.pressure_ratio is not None:
            mainvar.funclist.append(self.MomentumB)
        if self.power is not None:
            mainvar.funclist.append(self.PowerB)
        mainvar.streaminletlist.append(self.inlet)                              # Adding the inlet stream index into the mainvar.streaminletlist
        mainvar.streamoutletlist.append(self.outlet)                            # Adding the outlet stream index into the mainvar.streamoutletlist
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
        return
    def EnergyB(self):
        return mainvar.streamlist[self.outlet].h - mainvar.streamlist[self.inlet].h - (1/self.eta_is) * ( RF.H_PS(mainvar.streamlist[self.outlet].z, mainvar.streamlist[self.outlet].p, mainvar.streamlist[self.inlet].s) - mainvar.streamlist[self.inlet].h )        
    def MomentumB(self):
        return mainvar.streamlist[self.outlet].p / mainvar.streamlist[self.inlet].p - self.pressure_ratio 
    class funcinitialize:                                                       # Class assignement of CompositionB function is required to avoid passing of additional variables to a function in the mainvar.funclist
        def __init__(self, i, inlet, outlet):
            self.i = i
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.CompositionB)   
        def CompositionB(self):
            return mainvar.streamlist[self.inlet].z[self.i] - mainvar.streamlist[self.outlet].z[self.i]
    def MassB(self):
        return mainvar.streamlist[self.inlet].m - mainvar.streamlist[self.outlet].m
    def PowerB(self):
        return mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m - mainvar.streamlist[self.inlet].h * mainvar.streamlist[self.inlet].m - self.power
    def PowerF(self):
        return mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m - mainvar.streamlist[self.inlet].h * mainvar.streamlist[self.inlet].m
    def MomentumF(self):
        return mainvar.streamlist[self.outlet].p / mainvar.streamlist[self.inlet].p
    def ExergyF(self):
        return mainvar.streamlist[self.inlet].m * (mainvar.streamlist[self.inlet].ex - mainvar.streamlist[self.outlet].ex) + self.power
    def __update__(self):                                                       # Updating characteristics of component after simulation
        self.pressure_ratio = self.MomentumF()
        self.power = self.PowerF()
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return
    
    
class Turbine:
    """
    Input paramters:
        inlet_stream_number - index of inlet cycle stream starting with 1;
        outlet_stream_number - index of outlet cycle stream starting with 1;
        eta_is - isentropic efficiency (1 by default);
        pressure_ratio - p_inlet / p_outlet (optional);
        power - output power in kW (optional).
    Added equations:
        EnergyB - energy balance (h_in - h_out - (1/eta_is)*(h_in - h_outs))
        MassB - mass balance (m_in - m_out)        
        Optionally:
            CompositionB - composition equilibrium for each of mixture components except the last one (z_in[i] - z_out[i]), where i - component of the mixture; not defined for a pure fluid or if composition is fixed
            MomentumB - pressure balance if pressure ratio is given (p_in/p_out - pressure_ratio)
            PowerB - full energy balance if power is given (h_in*m_in - h_out*m_out - power)       
    PowerF, ExergyF, MomentumF - available functions to calculate actual power, exergy loss and pressure ratio
    """
    global mainvar
    def __init__(self, inlet_stream_number, outlet_stream_number, eta_is = 1, pressure_ratio = None, power = None):
        self.name = "Turbine"
        self.eta_is = eta_is
        self.inlet = inlet_stream_number - 1                                    # Internal index of inlet stream starting with 0
        self.outlet = outlet_stream_number - 1                                  # Internal index of outlet stream starting with 0
        self.pressure_ratio = pressure_ratio
        self.power = power
        self.exergy_loss = None
        mainvar.funclist.append(self.EnergyB)
        mainvar.funclist.append(self.MassB)  
        if mainvar.fixcomposition is None:
            for i in range(mainvar.fluid.componentsnumber-1):
                Turbine.funcinitialize(i, self.inlet, self.outlet)              # Adding the composition equilibrium functions for each of mixture components except the last one
        if self.pressure_ratio is not None:
            mainvar.funclist.append(self.MomentumB)
        if self.power is not None:
            mainvar.funclist.append(self.PowerB)
        mainvar.streaminletlist.append(self.inlet)                              # Adding the inlet stream index into the mainvar.streaminletlist
        mainvar.streamoutletlist.append(self.outlet)                            # Adding the outlet stream index into the mainvar.streamoutletlist 
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
        return
    def EnergyB(self):
        return mainvar.streamlist[self.inlet].h - mainvar.streamlist[self.outlet].h - self.eta_is * (mainvar.streamlist[self.inlet].h - RF.H_PS(mainvar.streamlist[self.outlet].z, mainvar.streamlist[self.outlet].p, mainvar.streamlist[self.inlet].s))        
    def MassB(self):
        return mainvar.streamlist[self.inlet].m - mainvar.streamlist[self.outlet].m
    class funcinitialize:                                                       # Class assignement of CompositionB function is required to avoid passing of additional variables to a function in the mainvar.funclist
        def __init__(self, i, inlet, outlet):
            self.i = i
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.CompositionB)   
        def CompositionB(self):
            return (mainvar.streamlist[self.inlet].z[self.i] - mainvar.streamlist[self.outlet].z[self.i])
    def MomentumB(self):
        return mainvar.streamlist[self.inlet].p / mainvar.streamlist[self.outlet].p - self.pressure_ratio 
    def PowerB(self):
        return mainvar.streamlist[self.inlet].h * mainvar.streamlist[self.inlet].m - mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m - self.power
    def PowerF(self):
        return mainvar.streamlist[self.inlet].h * mainvar.streamlist[self.inlet].m - mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m
    def MomentumF(self):
        return mainvar.streamlist[self.inlet].p / mainvar.streamlist[self.outlet].p
    def ExergyF(self):
        return mainvar.streamlist[self.inlet].m * (mainvar.streamlist[self.inlet].ex - mainvar.streamlist[self.outlet].ex) - self.power
    def __update__(self):                                                       # Updating characteristics of component after simulation
        self.pressure_ratio = self.MomentumF()
        self.power = self.PowerF()
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return


class Cell:  
    """
    Class to define one side of a multiflow heat exchanger.
    Each cell should have one inlet and one outlet. Heat exchangers with intermediate stream addition need to be splitted into simple heat exhcnagers.
    """
    def __init__(self, inlet, outlet, pressure_drop = 0, capacity = None):
        self.inlet = inlet - 1                                                  # Internal index of inlet stream starting with 0
        self.outlet = outlet - 1                                                # Internal index of outlet stream starting with 0
        self.pressure_drop = pressure_drop                                      # Pressure drop in bar (p_in - p_out)
        self.capacity = capacity                                                # Heat capacity of a heat exchanger side in kW
        return
           
class Heat_exchanger:
    """
    Multiflow heat exchanger with 2 or more sides
    Input parameters:
        cells - number of specified flows (sides) of a heat exchanger; shold be specified as a list of Cell classes: [Cell(...),Cell(...),]
        NTU - if required, the Number of Transfer Units can be spesified for a heat exhcnager. The outlet temperature of a cold stream will be set to meet the NTU specification (see function TC_NTU)
        cell1_for_NTU, cell2_for_NTU - NTU is calculated from inlet and outlet temperatures of heat exchanger streams, thus, for a heat exhcnager with more than 2 sides the cells for NTU calculations need to be specified (by default first two cells will be chosen).
    Added equations:
        EnergyB - full energy balance sum(h_in*m_in - h_out*m_out) for all cells
        MassB - mass balance for each cell (m_in - m_out)
        MomentumB - pressure balance for each cell (p_in - p_out)
        Optionally:
            CapacityB - power balance for each stream if capacity is known (h_in*m_in - h_out*m_out - capacity)
            CompositionB - composition equilibrium for each of mixture components except the last one (z_in[i] - z_out[i]), where i - component of the mixture; not defined for a pure fluid or if composition is fixed
            NTUB - equality of NTU to specified value (defined through the outlet temperature of a cold stream for a counterflow heat exchanger)
    ExergyF - available function to calculate exergy losses
    HXdef - function to calculate integral NTU value, pinch point and print Q-T diagram
    """
    global mainvar
    def __init__(self, cells: list, NTU = None, cell1_for_NTU = None, cell2_for_NTU = None):
        self.name = "Heat exchanger"
        self.number_of_cells = len(cells)                                       # Number of specified flows (sides) of multiflow heat exchanger
        self.cells = cells                                                      # List of heat exchanger sides (list of Cell classes)
        self.NTU = NTU                                                          # Desired Number of Transfer Units
        self.exergy_loss = None
        self.cell1_for_NTU = 0 if (cell1_for_NTU is None) else (cell1_for_NTU - 1)  # Setting the first cell number for NTU calculation
        self.cell2_for_NTU = 1 if (cell2_for_NTU is None) else (cell2_for_NTU - 1)  # Setting the second cell number for NTU calculation
        mainvar.funclist.append(self.EnergyB)   
        for i in range(self.number_of_cells):
            mainvar.streaminletlist.append(self.cells[i].inlet)                 # Adding the inlet stream index into the mainvar.streaminletlist
            mainvar.streamoutletlist.append(self.cells[i].outlet)               # Adding the outlet stream index into the mainvar.streamoutletlist
            Heat_exchanger.cellfuncinitialize(i,self.cells)
            if mainvar.fixcomposition is None:
                for j in range(mainvar.fluid.componentsnumber-1):
                    Heat_exchanger.zfuncinitialize(i,j,self.cells) 
        if self.NTU is not None:
            mainvar.funclist.append(self.NTUB)           
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
    def EnergyB(self):
        f = 0
        for i in range(self.number_of_cells):
            f = f + ( mainvar.streamlist[self.cells[i].inlet].h * mainvar.streamlist[self.cells[i].inlet].m - mainvar.streamlist[self.cells[i].outlet].h * mainvar.streamlist[self.cells[i].outlet].m )
        return f      
    class cellfuncinitialize:
        def __init__(self,i,cells):
            self.i = i
            self.cells = cells
            mainvar.funclist.append(self.MomentumB)
            mainvar.funclist.append(self.MassB)
            if self.cells[self.i].capacity is not None:
                mainvar.funclist.append(self.CapacityB) 
        def MassB(self): 
            return mainvar.streamlist[self.cells[self.i].inlet].m - mainvar.streamlist[self.cells[self.i].outlet].m
        def MomentumB(self):   
            return mainvar.streamlist[self.cells[self.i].inlet].p - mainvar.streamlist[self.cells[self.i].outlet].p - self.cells[self.i].pressure_drop
        def CapacityB(self):        
            return mainvar.streamlist[self.cells[self.i].inlet].m * mainvar.streamlist[self.cells[self.i].inlet].h - mainvar.streamlist[self.cells[self.i].outlet].m * mainvar.streamlist[self.cells[self.i].outlet].h - self.cells[self.i].capacity 
    class zfuncinitialize:
        def __init__(self,i,j,cells):
            self.i = i
            self.j = j
            self.cells = cells
            mainvar.funclist.append(self.CompositionB) 
        def CompositionB(self): 
            return (mainvar.streamlist[self.cells[self.i].inlet].z[self.j] - mainvar.streamlist[self.cells[self.i].outlet].z[self.j])
    def ExergyF(self):
        f = 0
        for i in range(self.number_of_cells):
            f = f + mainvar.streamlist[self.cells[i].inlet].m * ( mainvar.streamlist[self.cells[i].inlet].ex - mainvar.streamlist[self.cells[i].outlet].ex )
        return f 
    def NTUB(self):
        return TC_NTU(mainvar.streamlist[self.cells[self.cell1_for_NTU].inlet].T,mainvar.streamlist[self.cells[self.cell1_for_NTU].outlet].T,mainvar.streamlist[self.cells[self.cell2_for_NTU].inlet].T,mainvar.streamlist[self.cells[self.cell2_for_NTU].outlet].T,self.NTU) - mainvar.streamlist[self.cells[self.cell2_for_NTU].outlet].T      
    def HXdef(self, Para='mindT', Pval=None, savedata=False, plot=False):
        """
        Function to determine either the pinch point or an integral number of transfer units (NTU) of a multistream heat exchnager and its difference from a desired value (Pval).
        Additionally a Q-T diagram can be printed.
        """
        def HXplot(hxcurves):
            """
            Function to plot Q-T diagram
                hxcurves - total matrix of entalpies and temperature of cold and warm streams
            """
            import matplotlib.pyplot as plt
            y1 = hxcurves[:,1]
            y2 = hxcurves[:,2]
            x = hxcurves[:,0]
            fig = plt.figure()
            ax = fig.add_subplot(1, 1, 1)
            ax.plot(x,y1, color ='tab:blue')
            ax.plot(x,y2, color ='tab:orange')
            ax.set_ylabel("Temperature, K")
            ax.set_xlabel("Heat, kW")
            plt.show()
            return
        warmcells = []                                                          # List of all warm cells
        coldcells = []                                                          # List of all cold cells
        tempw = []                                                              # List of temperatures of warm side
        tempc = []                                                              # List of temperatures of cold side
        duty=0
        for i in range(self.number_of_cells):                                   # Identify warm and cold cells
            if mainvar.streamlist[self.cells[i].inlet].T < mainvar.streamlist[self.cells[i].outlet].T:
                coldcells.append(i)
                duty += mainvar.streamlist[self.cells[i].inlet].h * mainvar.streamlist[self.cells[i].inlet].m - mainvar.streamlist[self.cells[i].outlet].h * mainvar.streamlist[self.cells[i].outlet].m
                tempc.append(mainvar.streamlist[self.cells[i].inlet].T)
                tempc.append(mainvar.streamlist[self.cells[i].outlet].T)
            else:
                warmcells.append(i)
                tempw.append(mainvar.streamlist[self.cells[i].inlet].T)
                tempw.append(mainvar.streamlist[self.cells[i].outlet].T)
        #Fill the list of temperatures with intermediate temperatures
        tempc = np.unique(np.append(np.linspace(np.min(tempc), np.max(tempc), num=100), tempc))
        tempw = np.unique(np.append(np.linspace(np.min(tempw), np.max(tempc), num=100), tempw))
        hcsteps = np.zeros_like(tempc)                                          # List of enthalpies in the cold streams
        hwsteps = np.zeros_like(tempw)                                          # List of enthalpies in the warm streams
        for j in coldcells:                                                     # For each temperature step in the cold cells: sum up the enthalpy changes per temperature step
            z = mainvar.streamlist[self.cells[j].inlet].z
            hinlet = mainvar.streamlist[self.cells[j].inlet].h
            houtlet = mainvar.streamlist[self.cells[j].outlet].h
            massflow = mainvar.streamlist[self.cells[j].inlet].m
            for i in range(len(tempc)):
                if tempc[i] > mainvar.streamlist[self.cells[j].inlet].T and tempc[i] < mainvar.streamlist[self.cells[j].outlet].T:
                    hcsteps[i] += (RF.H_TP(z, tempc[i], mainvar.streamlist[self.cells[j].inlet].p) - hinlet)*massflow
                elif tempc[i] >= mainvar.streamlist[self.cells[j].outlet].T:
                    hcsteps[i] += (houtlet-hinlet) * massflow
        for j in warmcells:                                                     # For each temperature step in the warm cells: sum up the enthalpy changes per temperature step
            z = mainvar.streamlist[self.cells[j].inlet].z
            hinlet = mainvar.streamlist[self.cells[j].inlet].h
            houtlet = mainvar.streamlist[self.cells[j].outlet].h
            massflow = mainvar.streamlist[self.cells[j].inlet].m
            for i in range(len(tempw)):
                if tempw[i] > mainvar.streamlist[self.cells[j].outlet].T and tempw[i] < mainvar.streamlist[self.cells[j].inlet].T:
                    hwsteps[i] += (RF.H_TP(z, tempw[i], mainvar.streamlist[self.cells[j].inlet].p) - houtlet) * massflow
                elif tempw[i] >= mainvar.streamlist[self.cells[j].inlet].T:
                      hwsteps[i] +=-(houtlet-hinlet) * massflow
        if abs(max(hwsteps) - max(hcsteps)) > 1e-4:
            print('Enthalpies do not match, difference:', max(hwsteps) - max(hcsteps))
        enthalpies = np.concatenate((hwsteps, hcsteps))
        enthalpies = np.unique(np.sort(enthalpies))
        # The next block removes close enthalpy values that might cause problems for the evaluation
        enthalpiesN = []
        for i in range(len(enthalpies) - 1):
            if enthalpies[i+1] - enthalpies[i] > 1e-6 and enthalpies[i] > 0:
                enthalpiesN = np.append(enthalpiesN, enthalpies[i])
        if enthalpiesN[0]!=0:
            enthalpiesN = np.append([0], enthalpiesN)
        enthalpiesN = np.append(enthalpiesN, enthalpies[-1])
        enthalpies = enthalpiesN
        twarmcomp = np.interp(enthalpies, hwsteps, tempw)
        tcoldcomp = np.interp(enthalpies, hcsteps, tempc)
        hxcurves = np.column_stack((enthalpies, tcoldcomp, twarmcomp))
        if savedata == True:
            np.savetxt('hxcurve.txt', hxcurves)
        if plot == True:
            HXplot(hxcurves)
        if Para=='NTU':
            NTUs=[]
            for i in range(1, len(twarmcomp)):
                try:
                    NTUs = np.append(NTUs, NTUF(twarmcomp[i], twarmcomp[i-1], tcoldcomp[i-1], tcoldcomp[i], hxtype="general"))
                except:
                    print('NTU err')
                    raise
            value = sum(NTUs)
            if Pval is not None:
                res = value - Pval
            else:
                res=None
        else: # Para=='mindT':
            dt1 = np.array(twarmcomp - tcoldcomp)
            value = np.amin(dt1)
            if Pval is not None:
                res = value - Pval
            else:
                res = None
        return value, res
    def __update__(self):                                                       # Updating characteristics of component after simulation
        for i in range(self.number_of_cells):
            self.cells[i].capacity = mainvar.streamlist[self.cells[i].outlet].m * mainvar.streamlist[self.cells[i].outlet].h - mainvar.streamlist[self.cells[i].inlet].m * mainvar.streamlist[self.cells[i].inlet].h
            self.NTU = NTUF(mainvar.streamlist[self.cells[self.cell1_for_NTU].inlet].T, mainvar.streamlist[self.cells[self.cell1_for_NTU].outlet].T, mainvar.streamlist[self.cells[self.cell2_for_NTU].inlet].T, mainvar.streamlist[self.cells[self.cell2_for_NTU].outlet].T, hxtype = "counterflow")
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return
    
"""
Simple functions to calculate number of transfer units (NTU) of a heat exchnager
"""
def NTUF(Twin, Twout, Tcin, Tcout, hxtype = "counterflow"):
    if hxtype is "counterflow":
        C = (Twin - Twout) / (Tcout - Tcin)   
        if C < 1:
            E = (Tcout - Tcin) / (Twin - Tcin)
        else:
            C = 1 / C
            E = (Twin - Twout) / (Twin - Tcin)
        R = ((1 - C*E)/(1 - E))
        if R < 0 or R == 0:
            NTU = None
        else:
            NTU = (1 / (1 - C)) * log(R)  
    elif hxtype is "general":
        NTU = max((Twin - Twout), (Tcout - Tcin))
        NTU = NTU / (Twin - Tcout)
        if NTU < 0:
            NTU = 10000
            print('NTU<0')
    else:
        print("NTU model is not specified")
        NTU = None
    return NTU  

def TC_NTU(T1,T2,T3,T4,NTU): 
    """
    Outlet temperature of the cold side of the heat exchanger with given number of transfer units (for a counterflow heat exchanger only)    
    """
    C = (T1 - T2) / (T4 - T3)       
    if C<1:
        E = (1 - exp((C - 1) * NTU)) / (1 - C*exp((C - 1) * NTU))
        Ecold = E
    else:
        C = 1 / C
        E = (1 - exp((C - 1) * NTU))/(1 - C*exp((C - 1) * NTU))
        Ecold = E*C
    T4 = T3 + Ecold * (T1 - T3)
    return T4  
"""
......................
"""

       
class Simple_heat_exchanger:
    """
    Simple one-flow heat exchanger with given outlet temperature or heat load (in kW)
    Input paramters:
        inlet_stream_number - index of inlet cycle stream starting with 1;
        outlet_stream_number - index of outlet cycle stream starting with 1;
        pressure_frop - pressure drop in bar, p_inlet - p_outlet (optional);
        capacity - heat load in kW (>0 for cooler and <0 for heater) (optional);
        Tout - output temperature in K (optional).
    Added equations:
        MomentumB - pressure balance (p_in - p_out)
        MassB - mass balance (m_in - m_out)        
        Optionally:
            EnergyB - energy balance if heat load (capacity) is given (capacity - (h_out*m_out - h_in*m_in))
            CompositionB - composition equilibrium for each of mixture components except the last one (z_in[i] - z_out[i]), where i - component of the mixture; not defined for a pure fluid or if composition is fixed
            TemperatureB - temperature balance if outlet temperature is given (T_out - value)       
    PowerF, ExergyF, MomentumF - available functions to calculate actual power, exergy loss and pressure ratio
    """
    global mainvar
    def __init__(self,inlet_stream_number,outlet_stream_number, pressure_drop = 0, capacity = None, Tout = None):
        self.name = "Simple heat exchanger"
        self.inlet = inlet_stream_number - 1                                    # Internal index of inlet stream starting with 0
        self.outlet = outlet_stream_number - 1                                  # Internal index of outlet stream starting with 0
        self.pressure_drop = pressure_drop
        self.capacity = capacity 
        self.exergy_loss = None
        self.outlet_temperature = Tout
        mainvar.funclist.append(self.MomentumB) 
        mainvar.funclist.append(self.MassB)         
        if self.capacity is not None:
            mainvar.funclist.append(self.EnergyB)
        elif self.outlet_temperature is not None:
            mainvar.funclist.append(self.TemperatureB)
        mainvar.streaminletlist.append(self.inlet)                              # Adding the inlet stream index into the mainvar.streaminletlist
        mainvar.streamoutletlist.append(self.outlet)                            # Adding the outlet stream index into the mainvar.streamoutletlist
        if mainvar.fixcomposition is None:
            for i in range(mainvar.fluid.componentsnumber-1):
                Simple_heat_exchanger.funcinitialize(i, self.inlet, self.outlet) 
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
        return
    def EnergyB(self):
        return self.capacity - (mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m - mainvar.streamlist[self.inlet].h * mainvar.streamlist[self.inlet].m)
    def EnergyF(self):
        return (mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m - mainvar.streamlist[self.inlet].h * mainvar.streamlist[self.inlet].m)
    def MomentumB(self):        
        return mainvar.streamlist[self.inlet].p - mainvar.streamlist[self.outlet].p - self.pressure_drop
    def MassB(self):        
        return mainvar.streamlist[self.inlet].m - mainvar.streamlist[self.outlet].m         
    class funcinitialize:
        def __init__(self, i, inlet, outlet):
            self.i = i
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.CompositionB) 
            return
        def CompositionB(self):
            return (mainvar.streamlist[self.inlet].z[self.i] - mainvar.streamlist[self.outlet].z[self.i])    
    def TemperatureB(self):          
        return mainvar.streamlist[self.outlet].T - self.outlet_temperature
    def ExergyF(self):
        if (mainvar.streamlist[self.inlet].T - mainvar.exergyzerostate[0] )<0.0001 or (mainvar.streamlist[self.outlet].T - mainvar.exergyzerostate[0] )<0.0001:
            Tm = mainvar.exergyzerostate[0] 
        elif (np.abs( mainvar.streamlist[self.inlet].T - mainvar.streamlist[self.outlet].T ) < 0.0001):
            Tm = mainvar.streamlist[self.inlet].T
        else:
            Tm = (mainvar.streamlist[self.outlet].T - mainvar.streamlist[self.inlet].T ) / (np.log(mainvar.streamlist[self.outlet].T / mainvar.streamlist[self.inlet].T ))
        return mainvar.streamlist[self.inlet].m * (mainvar.streamlist[self.inlet].ex - mainvar.streamlist[self.outlet].ex) + self.capacity * (1 - ( mainvar.exergyzerostate[0] / Tm))
    def __update__(self):                                                       # Updating characteristics of component after simulation
        self.capacity = self.EnergyF()
        self.outlet_temperature = mainvar.streamlist[self.outlet].T
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return
        
class Separator:
    """
    Simple flow separator (keeps the flow composition and properties constant)
    Input paramters:
        inlet - index of inlet cycle stream starting with 1;
        *outlet - index of outlet cycle streams starting with 1 (so many as required);
    Added equations:
        MassB - mass balance (m_in - sum(m_out[i]))  
        TemperatureB - temperature equality of all outlets to the inlet (T_in - T_out[i])
        MomentumB - pressure equality of all outlets to the inlet (p_in - p_out[i])
        Optionally:
            CompositionB - composition equilibrium for each of mixture components except the last one (z_in[i] - z_out[i]), where i - component of the mixture; not defined for a pure fluid or if composition is fixed    
    ExergyF - available function to calculate exergy loss
    """
    global mainvar
    def __init__(self,inlet,*outlets):
        self.name = "Separator"
        self.inlet = inlet - 1                                                  # Internal index of inlet stream starting with 0
        self.outlet = []
        self.exergy_loss = None
        for value in outlets:
            self.outlet.append(value - 1)                                       # Internal indexes of outlet streams starting with 0
        mainvar.streaminletlist.append(self.inlet)                              # Adding the inlet stream index into the mainvar.streaminletlist
        mainvar.funclist.append(self.MassB)
        for i in range(len(self.outlet)):
            mainvar.streamoutletlist.append(self.outlet[i])                     # Adding the outlet stream indexes into the mainvar.streamoutletlist
            Separator.funcinitialize(self.inlet, self.outlet[i]) 
        if mainvar.fixcomposition is None:
            for i in range(mainvar.fluid.componentsnumber-1):
                for j in range(len(self.outlet)):              
                    Separator.zfuncinitialize(i, self.inlet, self.outlet[j])
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
    def MassB(self):
        f = mainvar.streamlist[self.inlet].m
        for i in range(len(self.outlet)):
            f = f - ( mainvar.streamlist[self.outlet[i]].m )
        return f 
    class funcinitialize:
        def __init__(self, inlet, outlet):
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.TemperatureB)
            mainvar.funclist.append(self.MomentumB)  
            return
        def TemperatureB(self):          
            return mainvar.streamlist[self.inlet].T - mainvar.streamlist[self.outlet].T
        def MomentumB(self):          
            return mainvar.streamlist[self.inlet].p - mainvar.streamlist[self.outlet].p
    class zfuncinitialize:
        def __init__(self, i, inlet, outlet):
            self.i = i
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.CompositionB)
            return
        def CompositionB(self): 
            return mainvar.streamlist[self.inlet].z[self.i] - mainvar.streamlist[self.outlet].z[self.i] 
    def ExergyF(self):
        f = mainvar.streamlist[self.inlet].m * mainvar.streamlist[self.inlet].ex
        for i in range(len(self.outlet)):
            f = f - ( mainvar.streamlist[self.outlet[i]].m * mainvar.streamlist[self.outlet[i]].ex )
        return f
    def __update__(self):                                                       # Updating characteristics of component after simulation
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return


class Mixer:
    """
    Flow mixer (defined for streams with equal pressure)
    Input paramters:
        outlet - index of outlet cycle stream starting with 1;
        *inlets - index of inlet cycle streams starting with 1 (so many as required);
    Added equations:
        EnergyB - full energy balance (h_out*m_out - sum(h_in*m_in))
        MassB - full mass balance (m_out - sum(m_in[i]))  
        MomentumB - pressure equality between all the inlets and outlet stream (p_out - p_in[i])
        Optionally:
            CompositionB - composition equilibrium for each of mixture components except the last one (calculated from compositions and mass flows of inlet streams); not defined for a pure fluid    
    ExergyF - available function to calculate exergy loss
    """
    global mainvar
    def __init__(self,outlet,*inlets):
        self.name = "Mixer"
        self.outlet = outlet - 1                                                # Internal index of outlet stream starting with 0
        self.inlet = [] 
        self.exergy_loss = None
        for value in inlets:
            self.inlet.append(value - 1)                                        # Internal indexes of inlet streams starting with 0
        mainvar.streamoutletlist.append(self.outlet)                            # Adding the outlet stream index into the mainvar.streamoutletlist
        mainvar.funclist.append(self.EnergyB)
        mainvar.funclist.append(self.MassB)
        if mainvar.fixcomposition is None:
            for i in range(mainvar.fluid.componentsnumber-1):
                Mixer.zfuncinitialize(i,self.outlet,self.inlet)
        for i in range(len(self.inlet)):
            mainvar.streaminletlist.append(self.inlet[i])                       # Adding the inlet streams index into the mainvar.streaminletlist
            Mixer.funcinitialize(self.outlet, self.inlet[i])
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
    def MassB(self):
        f = mainvar.streamlist[self.outlet].m
        for j in range(len(self.inlet)):
            f = f - ( mainvar.streamlist[self.inlet[j]].m )
        return f
    def EnergyB(self):
        f = mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m
        for i in range(len(self.inlet)):
            f = f - ( mainvar.streamlist[self.inlet[i]].h * mainvar.streamlist[self.inlet[i]].m )
        return f 
    class zfuncinitialize:
        def __init__(self, i, outlet, inlet):
            self.i = i
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.CompositionB)        
        def CompositionB(self):
            f = np.array([0]*mainvar.fluid.componentsnumber, dtype=float)
            for j in range(len(self.inlet)):
                Zmass=RF.MassComposition(mainvar.streamlist[self.inlet[j]].z)[:mainvar.fluid.componentsnumber]
                Zmass = [item*mainvar.streamlist[self.inlet[j]].m for item in Zmass]
                f = f + ( Zmass)           
            Z = RF.MoleComposition(f / mainvar.streamlist[self.outlet].m)
            return mainvar.streamlist[self.outlet].z[self.i] - Z[self.i]
    class funcinitialize:
        def __init__(self, outlet, inlet):
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.MomentumB)
        def MomentumB(self):          
            return mainvar.streamlist[self.inlet].p - mainvar.streamlist[self.outlet].p  
    def ExergyF(self):
        f = (-1) * mainvar.streamlist[self.outlet].m * mainvar.streamlist[self.outlet].ex
        for i in range(len(self.inlet)):
            f = f + ( mainvar.streamlist[self.inlet[i]].m * mainvar.streamlist[self.inlet[i]].ex )
        return f
    def __update__(self):                                                       # Updating characteristics of component after simulation
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return


class Valve:
    """
    Input paramters:
        inlet_stream_number - index of inlet cycle stream starting with 1;
        outlet_stream_number - index of outlet cycle stream starting with 1;
        pressure_frop - pressure drop in bar, p_inlet - p_outlet (optional);
        outlet_pressure - outlet pressure in bar (optional).
    Added equations:
        EnergyB - full energy balance (h_in*m_in - h_out*m_out)
        MassB - full mass balance (m_in - m_out)  
        MomentumB1 / Momentum B2 - pressure balance depending on specified parameter: for pressure_drop (p_in - p_out - pressure_drop) or for outlet_pressure (p_out - value)
        Optionally:
            CompositionB - composition equilibrium for each of mixture components except the last one (calculated from compositions and mass flows of inlet streams); not defined for a pure fluid or if composition is fixed 
    ExergyF - available function to calculate exergy loss
    """
    global mainvar
    def __init__(self, inlet_stream_number, outlet_stream_number, outlet_pressure = None, pressure_drop = None):
        self.name = "Valve"
        self.inlet = inlet_stream_number - 1                                    # Internal index of inlet stream starting with 0
        self.outlet = outlet_stream_number - 1                                  # Internal index of outlet stream starting with 0
        self.pressure_drop = pressure_drop
        self.outlet_pressure = outlet_pressure
        self.exergy_loss = None
        mainvar.funclist.append(self.EnergyB)
        mainvar.funclist.append(self.MassB)        
        if self.outlet_pressure is not None:
            mainvar.funclist.append(self.MomentumB1)
        if self.pressure_drop is not None:
            mainvar.funclist.append(self.MomentumB2)
        mainvar.streaminletlist.append(self.inlet)                              # Adding the inlet stream index into the mainvar.streaminletlist
        mainvar.streamoutletlist.append(self.outlet)                            # Adding the outlet stream index into the mainvar.streamoutletlist
        if mainvar.fixcomposition is None:
            for i in range(mainvar.fluid.componentsnumber-1):
                Valve.funcinitialize(i, self.inlet, self.outlet) 
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
        return
    def EnergyB(self):
        return (mainvar.streamlist[self.inlet].h * mainvar.streamlist[self.inlet].m - mainvar.streamlist[self.outlet].h * mainvar.streamlist[self.outlet].m)
    def MassB(self):
        return (mainvar.streamlist[self.inlet].m - mainvar.streamlist[self.outlet].m)    
    class funcinitialize:
        def __init__(self, i, inlet, outlet):
            self.i = i
            self.inlet = inlet
            self.outlet = outlet
            mainvar.funclist.append(self.CompositionB)           
        def CompositionB(self):
            return (mainvar.streamlist[self.inlet].z[self.i] - mainvar.streamlist[self.outlet].z[self.i])
    def MomentumB1(self):
        return mainvar.streamlist[self.outlet].p - self.outlet_pressure 
    def MomentumB2(self):
        return mainvar.streamlist[self.inlet].p - mainvar.streamlist[self.outlet].p - self.pressure_drop
    def ExergyF(self):
        return mainvar.streamlist[self.inlet].m * (mainvar.streamlist[self.inlet].ex - mainvar.streamlist[self.outlet].ex)
    def __update__(self):                                                       # Updating characteristics of component after simulation
        self.outlet_pressure = mainvar.streamlist[self.outlet].p
        self.pressure_drop = mainvar.streamlist[self.inlet].p - mainvar.streamlist[self.outlet].p
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return



class PhaseSeparator:
    """
    Input paramters:
        inlet - index of inlet cycle stream starting with 1;
        outlet_vapor - index of outlet vapor stream starting with 1;
        outlet_liquid - index of outlet liquid stream starting with 1.
    Added equations:
        MassBL, MassBV - definition of mass flows of vapor (MassBV) and liquid (MassBL) stream from the composition and mass flow of the inlet stream (through rerfpropfunc.py)   
        MomentumBvap, MomentumBliq - pressure equality of vapor and liquid streams to the inlet stream
        TemperatureBvap, TemperatureBliq - temperature equality of vapor and liquid streams to the inlet stream
        Optionally:
            CompositionBV, CompositionBL - composition equilibrium for each of mixture components except the last one (calculated from the composition and mass flow of the inlet stream through rerfpropfunc.py); not defined for a pure fluid or if composition is fixed 
    ExergyF - available function to calculate exergy loss
    """    
    global mainvar
    def __init__(self, inlet, outlet_vapor, outlet_liquid):
        self.name = "Phase Separator"
        self.inlet = inlet - 1                                                  # Internal index of inlet stream starting with 0
        self.outlet_vapor = outlet_vapor - 1                                    # Internal index of outlet vapor stream starting with 0
        self.outlet_liquid = outlet_liquid - 1                                  # Internal index of outlet liquid stream starting with 0
        self.exergy_loss = None
        mainvar.streaminletlist.append(self.inlet)                              # Adding the inlet stream index into the mainvar.streaminletlist
        mainvar.streamoutletlist.append(self.outlet_vapor)                      # Adding the outlet stream index into the mainvar.streamoutletlist
        mainvar.streamoutletlist.append(self.outlet_liquid)                     # Adding the outlet stream index into the mainvar.streamoutletlist
        self.mvap = None                                                        # Mass flow of vapor outlet
        self.mliq = None                                                        # Mass flow of liquid outlet
        self.vapcomposition = None                                              # Composition of vapor outlet
        self.liqcomposition = None                                              # Composition of liquid outlet
        mainvar.funclist.append(self.MassBV)
        mainvar.funclist.append(self.MassBL)
        mainvar.funclist.append(self.MomentumBvap)
        mainvar.funclist.append(self.MomentumBliq)
        mainvar.funclist.append(self.TemperatureBvap)
        mainvar.funclist.append(self.TemperatureBliq)
        if mainvar.fixcomposition is None:
            for i in range(mainvar.fluid.componentsnumber-1):
                PhaseSeparator.funcinitialize(i, self.outlet_vapor, self.outlet_liquid, self.inlet) 
        mainvar.componentslist.append(self)                                     # Adding the initialized class into the list of components (mainvar.componentslist)
        return
    def MassBV(self):
        result = RF.vaporstring(mainvar.streamlist[self.inlet].T, mainvar.streamlist[self.inlet].s, mainvar.streamlist[self.inlet].m, mainvar.streamlist[self.inlet].z)
        self.mvap = result[0]
        self.vapcomposition = result[1]
        return mainvar.streamlist[self.outlet_vapor].m - self.mvap
    def MassBL(self):
        result = RF.liquidstring(mainvar.streamlist[self.inlet].T, mainvar.streamlist[self.inlet].s, mainvar.streamlist[self.inlet].m, mainvar.streamlist[self.inlet].z)
        self.mliq = result[0]
        self.liqcomposition = result[1]
        return mainvar.streamlist[self.outlet_liquid].m - self.mliq
    def MomentumBliq(self):
        return mainvar.streamlist[self.inlet].p - mainvar.streamlist[self.outlet_liquid].p
    def MomentumBvap(self):
        return mainvar.streamlist[self.inlet].p - mainvar.streamlist[self.outlet_vapor].p  
    def TemperatureBliq(self):
        return mainvar.streamlist[self.inlet].T - mainvar.streamlist[self.outlet_liquid].T
    def TemperatureBvap(self):
        return mainvar.streamlist[self.inlet].T - mainvar.streamlist[self.outlet_vapor].T
    def ExergyF(self):
        return mainvar.streamlist[self.inlet].m * mainvar.streamlist[self.inlet].ex - mainvar.streamlist[self.outlet_vapor].m * mainvar.streamlist[self.outlet_vapor].ex - mainvar.streamlist[self.outlet_liquid].m * mainvar.streamlist[self.outlet_liquid].ex
    class funcinitialize:
        def __init__(self, i, outlet_vapor, outlet_liquid, inlet):
            self.outlet_vapor = outlet_vapor
            self.outlet_liquid = outlet_liquid 
            self.i = i
            self.inlet = inlet
            mainvar.funclist.append(self.CompositionBV)
            mainvar.funclist.append(self.CompositionBL)
            self.mvap = None
            self.mliq = None
            self.vapcomposition = None
            self.liqcomposition = None
            return
        def CompositionBV(self):
            result = RF.vaporstring(mainvar.streamlist[self.inlet].T, mainvar.streamlist[self.inlet].s, mainvar.streamlist[self.inlet].m, mainvar.streamlist[self.inlet].z)
            self.mvap = result[0]
            self.vapcomposition = result[1]
            return (mainvar.streamlist[self.outlet_vapor].z[self.i] - self.vapcomposition[self.i])
        def CompositionBL(self):
            result = RF.liquidstring(mainvar.streamlist[self.inlet].T, mainvar.streamlist[self.inlet].s, mainvar.streamlist[self.inlet].m, mainvar.streamlist[self.inlet].z)
            self.mliq = result[0]
            self.liqcomposition = result[1]
            return (mainvar.streamlist[self.outlet_liquid].z[self.i] - self.liqcomposition[self.i])
    def __update__(self):                                                       # Updating characteristics of component after simulation
        if mainvar.exergyflag == True:
            self.exergy_loss = self.ExergyF()
        return


    
"""
.........................................................................
Solution finding

To find the parameters of cycle, the mainfunction() should be solved. 
    - This function takes an array of current variables from the solver, separates them into temperature, mass flow, pressure/entropy and composition.
    - Then it defines (updates) parameters of all the streams is mainvar.streamlist.
    - Then it finds the solution of all the functions in mainvar.funclist for current parameters of streams in mainvar.streamlist and returns them as an array to the solver.
    - The solution is found when all the functions in mainvar.funclist are equal to 0 (with a tolerance specified by the solver).

To setup the solver and solve mainfunction() the solution() function is used.
    - The solution() function specifies the type of solution (temperature-entropy or temperature-pressure based)
    - Then it calls the setup.solversetup() to define the initial values (mainvar.guess) and check whether the number of defined system equations is right.
    - Then the scipy.optimize.root module is calles to find the solution of mainfunction() using mainvar.guess as initial values.
    - solution() returns the resulting variable array and the boolean flag showing whether the calculation is successful.
"""

def mainfunction(variable) :
    global mainvar
    Temperature = variable[0:(mainvar.streamnumber)]                            # Splitting the current solver variable array into physical stream parameters
    Massflow = variable[(2*mainvar.streamnumber):(3*mainvar.streamnumber)]
    if mainvar.fixcomposition is not None:
        Composition = mainvar.fixcomposition                                    # If known composition is defined by the user, it is taken from maincar.fixcomposition
    else:
        if mainvar.fluid.componentsnumber is 1:                                 # If pure fluid is used, the [1.] composition array is defined
            Composition = [[1.]*1]*mainvar.streamnumber
        else:
            Z = variable[(3*mainvar.streamnumber):]
            Z = np.split(Z, mainvar.streamnumber)
            Composition = []
            for i in range(mainvar.streamnumber):                               # For a mixture the full array of compositions is defined (mainvar.guess takes the mole fractions of all mixture components except the last one; the last mole fraction is calculated as (1-sum(Z[i])))
                Composition.append(np.concatenate((Z[i],[1-sum(Z[i])])))
    if mainvar.stype=="TP":                                                     # Parameters of each stream are defined/updated. Please note: the class stream is each time called as __init__(), so mainvar.streamlist always gets new memory adress of its streams.
        Pressure = variable[mainvar.streamnumber:(2*mainvar.streamnumber)]
        mainvar.streamlist = [stream(Temperature[i], Massflow[i], Composition[i], i, p=Pressure[i]) for i in range(mainvar.streamnumber)]
    else:
        Entropy=variable[mainvar.streamnumber:(2*mainvar.streamnumber)]
        mainvar.streamlist = [stream(Temperature[i], Massflow[i], Composition[i], i, s=Entropy[i]) for i in range(mainvar.streamnumber)]           
    Erg=[]
    for f in mainvar.funclist:                                                  # Finding solution of each equation of the defined system
        Erg=np.concatenate((Erg,f()),axis=None)
    return Erg


def solution(stype="TS"):
    """
    "TP" - temperature and pressure based calculation (faster for calculations in gas phase)
    "TS" - temperature and entropy based calculations (for 2-phase calculations and ideal gas behavior)
    """
    global mainvar
    mainvar.stype = stype
    if mainvar.setupflag is False:
        setup.solversetup()                                                     # Setup of mainvar.guess and check of function number
    solution = opt.root(mainfunction, mainvar.guess, method='hybr')             # Solving the system of equations with scipy.optimize.root
    ier = solution.success                                                      # Flag of solver success
    mesg = solution.message                                                     # Error message
    solution = solution.x                                                       # Solution as an array of T,P,M,*Z for "TP" symulation and T,S,M,*Z  for "TS" simulation
    if ier == 1:    
        mainvar.guess = solution                                                # Updating all the components for found solution
        mainfunction(solution)
        for component in mainvar.componentslist:
            component.__update__()
    else:
        print(mesg)
    return solution, ier


"""
.........................................................................
Output
"""

class output:
    global mainvar
    
    def txtconsole(solution):  
        """
        Full output of resulting component and stream parameters
            solution - full output of function solution()
        """
        solution, ier = solution                                                
        print("*************************************************************")    
        for component in mainvar.componentslist:                                # Printing all the properties of all specified cycle components
            attrs = vars(component)
            for item in attrs.items():
                if item[0]=="cells":
                    for obj in item[1]:
                        attrs2=vars(obj)
                        for item in attrs2.items():
                            print("%s: %s" % item)
                elif (item[0] is "PHS"):
                    pass
                else:
                    print("%s: %s" % item)
            print()        
        Temperature=solution[0:(mainvar.streamnumber)]                          # Printing the solution for Temperatures, Pressures, Massflows and composition
        Massflow = solution[2*mainvar.streamnumber:(3*mainvar.streamnumber)]
        if mainvar.fixcomposition is not None:
            Composition = mainvar.fixcomposition
        else:
            Composition = []
            Z = solution[(3*mainvar.streamnumber):]
            for i in range(mainvar.streamnumber):
                ZZ=Z[i*(mainvar.fluid.componentsnumber-1) : (i+1)*(mainvar.fluid.componentsnumber-1)]
                Composition.append(np.concatenate((ZZ,[1-sum(ZZ)])))
        if mainvar.stype=="TP":
            Pressure=solution[mainvar.streamnumber:(2*mainvar.streamnumber)]
        else:
            Entropy=solution[mainvar.streamnumber:(2*mainvar.streamnumber)]
            Pressure=np.array([None]*mainvar.streamnumber, dtype=np.float)
            for i in range(mainvar.streamnumber):
                Pressure[i]=RF.P_TS(Composition[i],Temperature[i],Entropy[i])
        print('Temperatures:')
        print(', '.join(str(item) for item in Temperature))
        print()
        print('Pressures:') 
        print(', '.join(str(item) for item in Pressure))
        print()
        print('Massflows:') 
        print(', '.join(str(item) for item in Massflow))
        print()
        print('Compositions:') 
        print(', '.join(str(item) for item in Composition))
        print("*************************************************************")
        print("Success:",ier)
        return
    
    def exergylosses():
        """
        Output of exergy balance in form: "component name : exergy loss in kW, ratio to total losses in %"
        """
        print("*************************************************************")
        print("Exergetic losses:")
        exloss=0
        for component in mainvar.componentslist:
            exloss=exloss+component.exergy_loss
        for component in mainvar.componentslist:
            print(component.name, ":", component.exergy_loss, ",", np.round(component.exergy_loss/exloss*100, 2),"%")
        print("*************************************************************")
        return
    
    def excel(filename="Default"):
        """
        Saving of an excel protocol into the user program location path
        """
        global mainvar  
        import xlsxwriter
        import datetime
        print('Writing to excel...')    
        if filename=="Default":
            filename=datetime.datetime.now().strftime("%Y%m%d_%H%M%S") + '_CryoSolver_Protocol'+'.xlsx'
        workbook = xlsxwriter.Workbook(filename)
        worksheet = workbook.add_worksheet("Simulation Protocol")
        cellsformat = workbook.add_format({'bg_color': '#99CC00','bold': True})                             
        worksheet.set_row(0, cell_format = cellsformat)
        worksheet.write(0, 0, "Stream properties") 
        worksheet.write(1, 0, "Index")
        worksheet.write(2, 0, "T, K")
        worksheet.write(3, 0, "P, bar")
        worksheet.write(4, 0, "S, kJ/(kg*K)")
        worksheet.write(5, 0, "H, kJ/kg")
        worksheet.write(6, 0, "Exergy, kJ/kg")
        worksheet.write(7, 0, "M, kg/s")
        worksheet.write(8, 0, "Fluid name")
        worksheet.write(9, 0, "Mole composition")
        row = 1
        column = 1
        for stream in mainvar.streamlist:
            worksheet.write(1, column, stream.index)
            worksheet.write(2, column, stream.T)
            worksheet.write(3, column, stream.p)
            worksheet.write(4, column, stream.s)
            worksheet.write(5, column, stream.h)
            worksheet.write(6, column, stream.ex)
            worksheet.write(7, column, stream.m)
            worksheet.write(8, column, mainvar.fluid.name)
            worksheet.write(9, column, str(stream.z))
            column+=1
        worksheet.set_row(11, cell_format = cellsformat)
        worksheet.write(11, 0, "Component properties") 
        row = 12
        column = 0
        for component in mainvar.componentslist:
            attrs = vars(component)
            for item in attrs.items():
                if item[0]=="cells":
                    for obj in item[1]:
                        attrs2=vars(obj)
                        for item in attrs2.items():
                            worksheet.write(row, column + 1, item[0])
                            if (type(item[1]) in (int, float, np.float64)):
                                worksheet.write(row, column + 2, item[1])
                            else:
                                worksheet.write(row, column + 2, str(item[1]))
                            row+=1
                elif (item[0] is "PHS"):
                    pass
                else:
                    worksheet.write(row, column, item[0])
                    if (type(item[1]) in (int, float, np.float64)):
                        worksheet.write(row, column + 2, item[1])
                    else:
                        worksheet.write(row, column + 2, str(item[1]))
                    row+=1
            row+=1
        workbook.close() 
        print("File saved:", filename)
        return
  

