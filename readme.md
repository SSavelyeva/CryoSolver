# CryoSolver. Package for cryogenic cycle simulation in Python
# Readme Ver. 1.0.0

## Introduction

This is the documentation for the _CryoSolver_ package for calculation and optimization of cryogenic cycles in Python. The _CryoSolver_ offers a collection of classes and functions to automatically determine and solve the system of equations such as mass and energy balances to define the parameters at cycle points. These equations are determined according to the cycle components and flows specified by the user. Its advantage is the possibility to add custom functions and components characteristics due to an open and flexible Python environment. The package does not offer a large variety of tools in contrast to the industrial software, but can be sufficient for research purposes.

Please note that the package is a by-product of the research activity and the authors are not responsible for any issues connected with its utilisation.

## Licensing

The _CryoSolver_ package is available under the [3-Clause BSD License](<https://opensource.org/licenses/BSD-3-Clause>) and can be used, copied, modified for free without permission on condition of an authorship acknowledgement:

> Copyright 2020 Sofiya Savelyeva, Steffen Klöppel (KKT, TU Dresden)

To cite the package, please reference to its [documentation](CryoSolver_documentation_Ver.1.0.0.pdf). 

## Installation

The package has been developed for _Python 3.6.4_ (or higher) and the [_Refprop 10.0_](https://www.nist.gov/srd/refprop). 

Note: It is possible to use CoolProp library instead of Refprop. The _refpropfunc.py_ needs to be respectively modified by the user. 

The [_ctRefprop_](https://pypi.org/project/ctREFPROP/) wrapper must be installed to allow the communication between Refprop and Python.

Other Python libraries used: _scipy_, _numpy_, _math_, _matplotlib_.

The simulations can be made with and without installation of the package.

### Import from a path:

The simulations can be made in an input file saved directly in the parent folder of the package (".../Cryosolver-1.0.0/").

Alternatively, to import the package in a random directory, the package path needs to be once added into the system path list:

```
import sys
sys.path.append("Path/to/CryoSolver-1.0.0/") 
```

### Global installation:

Before the global installation please execute any example file from the "/Cryosolver-1.0.0/examples/" folder and make sure that the program works properly. In case of errors with Refprop please refer to <https://www.nist.gov/srd/refprop> and modify the _refpropfunc.py_. In case of installation on OSX or Linux some changes in the Refprop setup can be required.

The installation can be made from the terminal using the package installer for Python ([pip](https://pip.pypa.io/en/stable/)): 

```
>pip install C:/.../CryoSolver-1.0.0
```

After that the package can be imported in any Python file with the command:

```
import cryosolver
```

## Simulation principle

The simulation principle and detailed package content is described in the [CryoSolver documentation file](CryoSolver_documentation_Ver.1.0.0.pdf).

_CryoSolver_ consists of two main files: 

* _refpropfunc.py_
* _solver.py_

The _refpropfunc.py_ contains simplified functions to calculate the fluid properties using the Refprop 10 library. The _solver.py_ is the main program containing the library of functions to determine cycle components, to setup parameters of flows and initial values of calculations and to solve the system of energy, mass balance and other specified equations. The composition Z should be given on a **mole basis** by default.

Units of measure used in the program:

|Name | Designation in code | units |
|-|-|-|
|Temperature | T | K | 
|Pressure | P | bar (abs) |  
|Enthalpy | H | kJ/kg | 
|Entropy | S | kJ/(kgK)| 
|Density | D | kg/m^3 |  
|Isobaric and isochoric heat capacity | Cp, Cv | kJ/(kgK) |	
|Mass flow | m | kg/s |	
|Power | power | kW |  
|Heat load | capacity | kW |
|Isentropic efficiency | eta_is | \- |

The composition Z should be given on a **mole basis** by default.

## Using the package

All the calculations should be preferably done in a separate input program file. Examples of such files are located in the /CryoSolver-1.0.0/examples/ folder.

### Main parts of the user program

The calculation is divided into the following main steps:

1. Solver import:

	```
	import cryosolver.solver as CS
	```

2. Reset of global variables.

	As global variable (_mainvar_ class object) is used by the solver, it needs to be created or reset before each new simulation to delete old values:

	```
	CS.start()
	```

	If necessary, by means of this function it is possible to change a preferable zero state for the exergy calculation. By default the zero state is at 300 K and 1 bar.

	```
	CS.start([303,1.013])
	```

	Furthermore, if a cycle with known mixture composition for all streams should be simulated, it is possible to put in here an array of compositions  to skip the calculation of composition balances (_x1_ and _x2_ - value of mole fractions of mixture components for each of 3 streams):


	```
	CS.start(fixcomposition=[[x1,x2],[x1,x2],[x1,x2]])
	```


3. Fluid specification:

	```
	CS.setup.fluidsetup("Nitrogen.FLD", "Oxygen.FLD")
	```

	Please refer to the [_CryoSolver_ documentation file](CryoSolver_documentation_Ver.1.0.0.pdf) and notes in the program code.

4. Assignement of cycle components by the user with class objects included into the package. 

	Depending on given component specification the required balance functions are added by the solver into the list of equations (_cryosolver.solver.mainvar.funclist_).

	Thus, for example, if the turbine with 85% efficiency is assigned between the 1st and 2nd streams,

	```
	CS.Turbine(1,2,0.85)
	```

	the following functions will be added to _cryosolver.solver.mainvar.funclist_:

	* Energy balance: <img src="https://render.githubusercontent.com/render/math?math=EnergyB()=h_{inlet}-h_{outlet}-\eta_s\cdot(h_{inlet}-h_{outlet} \big|_s)">

	* Mass balance: <img src="https://render.githubusercontent.com/render/math?math=MassB()=\dot{m}_{inlet}-\dot{m}_{outlet}">

	Here inlet and outlet indexes state for streams 1 and 2 respectively.

	In case if turbine power (here 700 kW) will be given:

	```
	CS.Turbine(1,2,0.85,power=700)
	```

	additional function will be added:

	* Energy balance: <img src="https://render.githubusercontent.com/render/math?math=PowerB()=h_{inlet}\cdot \dot{m}_{inlet} - h_{outlet}\cdot \dot{m}_{outlet} - power">

5. Specification of fluid composition. 

	This step should be done only for mixtures and if no fixed composition array is specified by _CS.start()_:

	```
	CS.setup.parameterset(1, "Z",[0.78,0.22]) 
	```

	6. Fixing the known parameters.

	The respective functions are added into the list. For example, to fix the temperature and pressure at the turbine outlet as 100 K, 1 bar, it should be written:

	```
	CS.setup.parameterset(2,"T",100)
	CS.setup.parameterset(2,"P",1)
	```

7. Additional functions can be added manually into the system. 

	For example, fixing the temperature difference through the turbine as 100 K: 

	```
	def userfunction():
		return CS.mainvar.streamlist[0].T - CS.mainvar.streamlist[1].T - 100	
	CS.mainvar.funclist.append(userfunction)
	```

	Please note, here the indexes are related to the internal solver array, so for the stream _i_ the index _i-1_ should be used (see rules below).	

8. Specification of an array of initial values.

	It depends on specific problem, how close they should be to the solution. Usually, for a simple case, any physically meaningful values could be assigned. The array should mandatory have 3 list elements for temperature, pressure and mass flow and an array of mole compositions except the last component of the mixture:

	```
	CS.setup.initvaluesset( ([300,200], [5,1], [1,1],[[0.78],[0.78]])  )
	```

	Here 0.78 is the initial mole fraction of nitrogen, the one of oxygen will be calculated as _(1 - 0.78)_.

	In case of a pure fluid or if the composition is fixed at the beginning, the initial values only for temperatures, pressures and mass flows need to be specified:

	```
	CS.setup.initvaluesset( ([300,200], [5,1], [1,1])  )
	```

9. Finding the solution and result output.

	For definition of an optimization/case study functions please refer to example files.

	```
	s = CS.solution("TP")
	CS.txtoutput(s)
	```

	It is possible to perform a simulation on temperature/pressure ("TP") or on temperature/entropy ("TS") based calculation. The second one is more stable for ideal gases and 2-phase region simulations. However, the initial values should be anyway setup as temperature/pressure based and will be automatically recalculated into the entropy. Please note, that manual printing of the solution in case of "TS" simulation will return the entropy values in the middle. In the _CS.textoutput()_ function they are already transformed back to show pressures.

	
The execution of the above code for nitrogen and oxygen (the full code is saved in /CryoSolver-1.0.0/examples/Example_0.py) will give the output:

```
Number of specified streams:  2
Specified open streams:                                        
inlets: [0] outlets: [1]
*************************************************************
name: Turbine
eta_is: 0.85
inlet: 0
outlet: 1
pressure_ratio: 15.112917789718212
power: 699.9999999999162
exergy_loss: 406.64579563013615

Temperatures:
200.0, 100.0

Pressures:
15.112917789718212, 1.0

Massflows:
5.093246893050264, 5.093246893050264

Compositions:
[0.6 0.4], [0.6 0.4]
*************************************************************
Success: True
```

By the specified open streams all the inlets and outlets of the cycle are shown regardless of whether they belong to a closed loop. This output allows to find possible errors in stream index definitions.


### Important simulation rules


#### Rule to define a closed cycle


As well as in most of industrial simulation software, it is not possible to directly setup a fully closed cycle. This would cause an over-definition of the equation system due to looping of the mass flow and composition balance from one stream to another. Thus, the cycle should be defined with at least **one opening** and two **closing equation of pressure and temperature balance**.

For comfort, the additional command including the setup of the pressure equality of two streams is added to the solver:

```
CS.setup.closeloop(1,6)
```
	
The temperature balance is not included, so temperatures can be specified at both ends if required.


#### Stream index definition

The streams need to be **numbered starting with 1** and not to have breaks in between. This is important for the solver to find the total number of variables.

However, the solver shifts these indexes to start the assignment of internal array elements with **0**.  Thus, to access the stream properties from the solver (for example, defining user functions with _cryosolver.solver.mainvar.streamlist_) **the (i-1) index** should be used for each stream  _i_.

Example: to set the temperature of the first cycle stream as 300 K using the module function:

```
CS.setup.parameterset(1,"T",300)
```

To print the temperature of the first stream from the _streamlist_ array:

```
print(CS.mainvar.streamlist[0].T)
```
	
## Solution failures

The solution will only converge if the problem statement allows to find a definitive solution. Thus, all the components and parameters need to be specified in such a way that they don't contradict with physical laws and specified functions governing all these variables. Moreover, the specification of initial values for a solver can have a high impact on the convergence. 

If no solution is found:
* check indexes of all specified streams and parameters (according to the stream assignment rule);
* check physical statement of the problem (for example, if solver finds a solution with negative values of variables, this most likely indicates incompatibility of set values);
* change initial values;
* for a more complex problem solver type and settings can be changed inside the _solver.py_ module;
* balance equations could be multiplied by some weight coefficients to change the tolerance of solution;
* variables could be multiplied by some scaling factors before being passed to the solver and then scaled back inside the _mainfunction()_ to change the solver sensitivity;
* in case of problems with fluid properties refer to the Refprop documentation <https://www.nist.gov/srd/refprop>. 

Additionally, there can be failures in the calculation due to insufficient accuracy of solver solution or optimisation of a function with many peaks. These issues need to be solved separately for each specific problem.

## Examples

The explanation of example files is given in the [CryoSolver documentation file](CryoSolver_documentation_Ver.1.0.0.pdf).

## Acknowledgement

This package was developed at the TU Dresden within the Work Package 4 of the EASITrain Project. EASITrain – European Advanced Superconductivity Innovation and Training. This Marie Sklodowska-Curie Action (MSCA) Innovative Training Networks (ITN) receives funding from the European Union’s H2020 Framework Programme under grant agreement no. 764879, <http://easitrain.web.cern.ch/>.

## References

1. Lemmon, E.W., Bell, I.H., Huber, M.L., McLinden, M.O. NIST Standard Reference Database 23: Reference Fluid Thermodynamic and Transport Properties-REFPROP, Version 10.0, National Institute of Standards and Technology, Standard Reference Data Program, Gaithersburg, 2018. <http://dx.doi.org/10.18434/T4JS3C>

2. Savelyeva, S., Kloeppel, S, Haberstroh, Ch. CryoSolver. Package for cryogenic cycle simulation in Python. Documentation Ver. 1.0.0. Dresden, 2018. <https://doi.org/10.5281/zenodo.4001668> 
