# -*- coding: utf-8 -*-
"""
Created on Wed Mar 11 09:25:46 2020

CryoSolver package
@author: Sofiya Savelyeva, Steffen Klöppel
TU Dresden 2020
"""
   
"""
.........................................................................

refprop.py
Simplified functions to calculate fluid properties with Refprop 10 functions (High-level API)
Information about the Refprop:
https://www.nist.gov/system/files/documents/2018/05/23/refprop10a.pdf

The following units are used:
    - T - temperature in K;
    - P - pressure in bar;
    - H - enthalpy in kJ/kg;
    - S - entropy in kJ/(kg*K);
    - D - density in kg/m^3;
    - Cp, Cv - isobaric and isochoric heat capacity in kJ/(kg*K);
    - m - mass flow, kg/s.
Composition Z should be indicated in mole units.

The fluid setup is performed only once using the SETUPdll routine of the Refprop for all the pure fluids of the cycle. 
The array of mole composition Z has the lengh of number of pure fluids indicated by setup (fluid.componentsnumber). For a pure fluid stream the rest mole fractions should be simply indicated as 0.

The reference state is set to 300 K and 1 bar (abs) substracting the corresponding reference value from the calculated absolute enthalpy and entropy.
The internal Refprop subroutine is not used to avoid deviations during the mixture calculations.

"""
import os
from ctREFPROP.ctREFPROP import REFPROPFunctionLibrary
RP = REFPROPFunctionLibrary(os.environ['RPPREFIX']+'\REFPROP.DLL')
RP.SETPATHdll(os.environ['RPPREFIX'])


class fluid:
    """
    This class determines the fluid which properties are calculated.
    
    name: 'name.FLD' for pure fluids or ('name1.FLD','name2.FLD') for mixtures (Windows)  
        note: for MAC the name is ‘:fluids:name.fld’, for Unix: ‘[full_path]/fluids/name.fld’
    componentsnumber: number of mixture components
    """
    def __init__(self, fluids : list):
        self.componentsnumber = len(fluids)                                     # Number of defined fluids
        if len(fluids)>1:
            self.name = '|'.join(fluids)                                        # Fluid names separation in form "name1.FLD|name2.FLD"
        else:
            self.name = fluids[0]
        r = RP.SETUPdll(self.componentsnumber, self.name, "HMX.BNC", "DEF")     # Fluid setup for Refprop
        if r.ierr != 0:                                                         # Fluid setup error check
                print(r.herr)
        return


"""
Fluid properties
"""

def MoleComposition(MassComposition):
    """
    Conversion of the mass flow of components into the mole composition
        MassComposition - array of mass fractions of each pure fluid 
        MoleComposition - array of mole fractions
    """
    MoleComposition = RP.XMOLEdll(MassComposition).xmol
    return MoleComposition


def MassComposition(MoleComposition):
    """
    Conversion of the mole composition into the mass composition
        MoleComposition - array of mole fractions
        MassComposition - array of mass fractions
    """
    MassComposition = RP.XMASSdll(MoleComposition).xkg      
    return MassComposition

def H_TP(Z, T, P):
    """
    Enthalpy defined by temperature and pressure (reference state at 300 K and 1 bar)
        Z - array of molar composition
        T, P - temperature and pressure    
        Units are specified above
    """
    H = RP.ABFLSHdll('TP', T, P*100, Z, 2).h - RP.ABFLSHdll('TP', 300, 100, Z, 2).h  
    return H

def S_TP(Z, T, P): 
    """
    Entropy defined by temperature and pressure (reference state at 300 K and 1 bar)
        Z - array of molar composition
        T, P - temperature and pressure    
        Units are specified above
    """
    S = RP.ABFLSHdll('TP', T, P*100, Z, 2).s - RP.ABFLSHdll('TP', 300, 100, Z, 2).s 
    return S

def H_PS(Z, P, S):
    """
    Enthalpy defined by pressure and entropy (reference state at 300 K and 1 bar)
        Z - array of molar composition
        P, S - pressure and entropy (reference state at 300 K and 1 bar)
        Units are specified above
    """ 
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    HH = RP.ABFLSHdll('PS', P*100, SS, Z, 0).h - RP.ABFLSHdll('TP', 300, 100, Z, 0).h 
    H = HH / MOL  
    return H

def P_TS(Z, T, S): 
    """
    Pressure defined by temperature and entropy
        Z - array of molar composition
        T, S - temperature and entropy (reference state at 300 K and 1 bar)
        Units are specified above
    """
    MOL = RP.WMOLdll(Z) # Molar mass
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s   
    PP = RP.ABFLSHdll('TS', T, SS, Z, 0).P 
    P = PP / 100
    return P

def H_TS(Z, T, S): 
    """
    Enthalpy defined by temperature and entropy (reference state at 300 K and 1 bar)
        Z - array of molar composition
        T, S - temperature and entropy (reference state at 300 K and 1 bar)
        Units are specified above
    """       
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    HH = RP.ABFLSHdll('TS', T, SS, Z, 0).h - RP.ABFLSHdll('TP', 300, 100, Z, 0).h 
    H = HH / MOL  
    return H

def D_TP(Z, T, P):
    """
    Density defined by temperature and pressure
        Z - array of molar composition
        T, P - temperature and pressure
        Units are specified above
    """  
    D = RP.ABFLSHdll('TP', T, P*100, Z, 2).D      
    return D

def T_PS(Z, P, S): 
    """
    Temperature defined by pressure and entropy
        Z - array of molar composition
        P, S - pressure and entropy (reference state at 300 K and 1 bar)
        Units are specified above
    """      
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    SS = S*MOL + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    T = RP.ABFLSHdll('PS', P*100, SS, Z, 0).T      
    return T

def D_PS(Z, P, S):
    """
    Density defined by pressure and entropy
        Z - array of molar composition
        P, S - pressure and entropy (reference state at 300 K and 1 bar)
        Units are specified above
    """     
    MOL=RP.WMOLdll(Z)                                                           # Molar mass
    SS = S*MOL + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    DD = RP.ABFLSHdll('PS', P*100, SS, Z, 0).D 
    D = DD*MOL
    return D

def Cv_TP(Z, T, P):
    """
    Isochoric heat capacity defined by temperature and pressure
        Z - array of molar composition
        T, P - temperature and pressure
        Units are specified above
    """  
    Cv = RP.ABFLSHdll('TP', T, P*100, Z, 2).Cv      
    return Cv

def Cp_TP(Z, T, P): 
    """
    Isobaric heat capacity defined by temperature and pressure
        Z - array of molar composition
        T, P - temperature and pressure
        Units are specified above
    """      
    Cp = RP.ABFLSHdll('TP', T, P*100, Z, 2).Cp       
    return Cp

def T_PH(Z, P, H): 
    """
    Temperature defined by pressure and enthalpy
        Z - array of molar composition
        P, H - pressure and enthalpy (reference state at 300 K and 1 bar)
        Units are specified above
    """       
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    HH = H*MOL + RP.ABFLSHdll('TP', 300, 100, Z, 0).h 
    T = RP.ABFLSHdll('PH', P*100, HH, Z, 0).T   
    return T

def D_PH(Z, P, H):
    """
    Density defined by pressure and enthalpy
        Z - array of molar composition
        P, H - pressure and enthalpy (reference state at 300 K and 1 bar)
        Units are specified above
    """      
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    HH = H*MOL + RP.ABFLSHdll('TP', 300, 100, Z, 0).h 
    DD = RP.ABFLSHdll('PH', P*100, HH, Z, 0).D
    D = DD*MOL
    return D

def P_TH(Z, T, H): 
    """
    Pressure defined by temperature and enthalpy
        Z - array of molar composition
        T, H - temperature and enthalpy (reference state at 300 K and 1 bar)
        Units are specified above
    """      
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    HH = H*MOL + RP.ABFLSHdll('TP', 300, 100, Z, 0).h 
    PP = RP.ABFLSHdll('TH', T, HH, Z, 0).P
    P = PP / 100
    return P

def P_HS(Z, H, S): 
    """
    Pressure defined by enthalpy and entropy
        Z - array of molar composition
        H, S - enthalpy and entropy (reference state at 300 K and 1 bar)
        Units are specified above
    """         
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    HH = H*MOL + RP.ABFLSHdll('TP', 300, 100, Z, 0).h 
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s   
    PP = RP.ABFLSHdll('HS', HH, SS, Z, 0).P 
    P = PP / 100
    return P

def S_PH(Z, P, H):
    """
    Entropy defined by pressure and enthalpy (reference state at 300 K and 1 bar)
        Z - array of molar composition
        P, H - pressure and enthalpy (reference state at 300 K and 1 bar)
        Units are specified above
    """       
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    HH = MOL*H + RP.ABFLSHdll('TP', 300, 100, Z, 0).h        
    SS = RP.ABFLSHdll('PH', P*100, HH, Z, 0).s - RP.ABFLSHdll('TP', 300, 100, Z, 0).s  
    S = SS / MOL      
    return S

"""
Stream properties in the two-phase region
"""

def xliq_TP(Z, T, P): 
    """
    Liquid phase mole composition defined by temperature and pressure
        Z - array of molar composition
        T, P - temperature and pressure
        Units are specified above
        xliq - array of mole fractions corresponding to the number of fluid components 
    """
    xliq = RP.ABFLSHdll('TP', T, P*100, Z, 2).x     
    return xliq

def xvap_TP(Z, T, P):   
    """
    Vapor phase mole composition defined by temperature and pressure
        Z - array of molar composition
        T, P - temperature and pressure
        Units are specified above
        xvap - array of mole fractions corresponding to the number of fluid components 
    """
    xvap = RP.ABFLSHdll('TP', T, P*100, Z,2).y     
    return xvap

def xliq_TS(Z, T, S): 
    """
    Liquid phase mole composition defined by temperature and entropy 
        Z - array of molar composition
        T, S - temperature and entropy (reference state at 300 K and 1 bar)
        Units are specified above
        xliq - array of mole fractions corresponding to the number of fluid components 
    """
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    xliq = RP.ABFLSHdll('TS', T, SS, Z, 0).x     
    return xliq

def xvap_TS(Z, T, S):   
    """
    Vapor phase mole composition defined by temperature and entropy
        Z - array of molar composition
        T, S - temperature and entropy (reference state at 300 K and 1 bar)
        Units are specified above
        xvap - array of mole fractions corresponding to the number of fluid components 
    """
    MOL = RP.WMOLdll(Z)                                                         # Molar mass
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    xvap = RP.ABFLSHdll('TS', T, SS, Z, 0).y     
    return xvap

def vaporstring(T, S, M, Z): 
    """
    Mass flow and mole composition of the vapor phase defined by temperature, entropy and total mass flow
        Z - array of molar composition
        T, S - temperature and entropy (reference state at 300 K and 1 bar)
        m - total stream mass flow, kg/s
        Units are specified above
    Returns:
        mflash - total mass flow of the vapor phase, kg/s
        flash - composition of the vapor phase on a mole basis
    """   
    MOL = RP.WMOLdll(Z)
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    Allprops = RP.ABFLSHdll('TS', T, SS, Z, 0)                                  # All properties of the stream
    Q = Allprops.q                                                              # Vapor quality on a molar basis
    if Q<0: Q = 0
    elif Q>1: Q = 1
    flash = Allprops.y                                                          # Composition of the vapor phase - array of mole fractions 
    MOLvap = RP.WMOLdll(flash)                                                  # Molar mass of the vapor phase
    mflash = M * (Q * MOLvap / MOL)                                             # Mass flows of the vapor phase
    return mflash, flash

def liquidstring(T, S, M, Z):
    """
    Mass flow and mole composition of the liquid phase defined by temperature, entropy and total mass flow
        T, S - temperature and entropy (reference state at 300 K and 1 bar)
        m - total stream mass flow, kg/s
        Units are specified above
    Returns:
        mflash - total mass flow of the liquid phase, kg/s
        flashx - composition of the liquid phase on a mole basis
    """     
    MOL = RP.WMOLdll(Z)                                                         # Molar mass 
    SS = MOL*S + RP.ABFLSHdll('TP', 300, 100, Z, 0).s
    Allprops = RP.ABFLSHdll('TS', T, SS, Z, 0)                                  # All properties of the stream
    Q = Allprops.q                                                              # Vapor quality on a molar basis
    if Q<0: Q = 0
    elif Q>1: Q = 1
    flashx = Allprops.x                                                         # Composition of the liquid phase - array of mole fractions
    flashy = Allprops.y                                                         # Composition of the vapor phase - array of mole fractions
    MOLvap = RP.WMOLdll(flashy)                                                 # Molar mass of the vapor phase
    mflash = M * (1 - Q * MOLvap / MOL)                                         # Mass flows of the vapor phase
    return mflash, flashx


def S_TQ(Z, T, Q):
    """
    Entropy defined by pressure and enthalpy (reference state at 300 K and 1 bar)
        Z - array of molar composition
        T, Q - temperature and vapor quality
        Units are specified above
    """       
    MOL = RP.WMOLdll(Z)                                                         # Molar mass     
    SS = RP.ABFLSHdll('TQ', T, Q, Z, 0).s - RP.ABFLSHdll('TP', 300, 100, Z, 0).s  
    S = SS / MOL      
    return S
 
