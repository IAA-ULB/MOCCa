

# Tantalus & Hephaestos
> Copyright W. Ryssens, P.H. Heenen & M. Bender
>Wednesday, 09. December 2020 11:58AM 

----------
### An introduction to Greek mythology

Or, perhaps more appropriate
> I have been given access to this repository, and have no clue where to start?!

__Tantalus__ is a mean-field code, representing the single-particle wavefunctions on a Lagrange mesh, written in FORTRAN.
__Hephaestos__ is a (collection of) Python modules/scripts that partially writes (parts of) the FORTRAN source code for Tantalus. 


Hephaestos, is in fact a code generator or preprocessor for Tantalus. For different choices regarding types of functionals and symmetries, Hephaestos generates a version of Tantalus that is able to perform calculations for you, subject to the input you provided Hephaestos.

-------
### File structure

**The chief files** for compilation are

*   `Hephaestos.py`: Main python script of Hephaestos.
*  `Makefile`: Makefile to compile Tantalus. 

**Source code **folders are

*  `src_heph/` :  contains modules for Hephaestos (in python). 
*  `src_orig/` :  contains source code files (in FORTRAN) of Tantalus, to be completed by Hephaestos.
*  `src/` :   directory of the compileable Tantalus source code, where Hephaestos places the processed .f90 files.

Information on the various **functionals** are

*  `functionals/` : Contains the definition of various forms of functionals (`.func` files), as input for Hephaestos.
*  `parameterizations/` : Contains various parameterization file (`.param` files), each linked to a particular `.func` file. 

**Auxiliary folders** are

* `configs/` : Predefined sets of compilation instructions can be found here. 
* `exec/` :  The Makefile places the compiled executables here.
* `mod/` : Auxiliary folder where the Makefile places the compiled `.mod` files.
* `obj/` :  Auxiliary folder fo the compiled `.obj` files. 

**Other**

* `tantalus_examples/` : contains several examples on how to use the code, once compiled correctly. These examples can also be used as unit tests; reference output and a comparison script are provided.

-------
### How do I compile the code?

The easiest way is to type

	make

or 

	make CONFIG=default

1.  Runs Hephaestos with input defined in `configs/default.py`. 
      This particular file dictates:
	* a pre-defined type of functional, in this case `NLO.func`.
	* a symmetry choice, in this case the EV8 = maximally symmetric one.
2.  Final processed Fortran files end up in `src/`.
3.  Said files are compiled.
4. The resulting exe (Tantalus.default.exe in this case) is placed exe into `exec/`.


Currently defined options to pass into the Makefile are

* `CXX`: your preferred choice of compiler. Default = `gfortran`.
* `CXXFLAGS`:  options for your preferred compiler.  
* `CONFIG`:  name of the configuration file dictating the options to Hephaestos. 
      Note: 
	* without the .py file extension
	* "." cannot be part of the filename 

__Things do check if the code does not want to compile__

* Are the `src/`,`mod/` and `obj/` directories present in your directory tree?

-------
### How do I run Hephaestos?

Answer:	

Provided you have **python-3** installed
		
		python Hephaestos.py configs/default.py

where `configs/default.py` can be replaced with another configuration file.  Specifying a configuration file is __not optional__! The code will process files from the `src_orig/ ` directory and place processed files in the `src` directory. 

Separate running of Hephaestos is not immediately useful, as after that Tantalus still needs to be compiled. 

-------
### How do I run Tantalus? 

Answer

	./Tantalus.default.exe < tant.data

which will produce output to STDOUT and a .wf file at the end of the calculation. The contents of tant.data contain all the relevant input. A minimal example for such file is

	cat << EOF > tant.data
	# Definition of the nucleus 
	&nucleus
	neutrons=8, protons=8
	/
	# Parameters of the Lagrange mesh. 
	&mesh
	nx=12, ny=12, nz=12, dx=1.0
	/
	# The code will look, on a file forces.param, for the parameterization with 
	# this name.
	&func
	name_param='SLy4'
	/
	# Options for the pairing.
	&pairing
	/
	# maxiter = Maximum number of iterations to be performed
	&evolution
	maxiter=100
	/
	&scfiteration
	/
	# Number of neutron (nwn) and proton (nwp) spwfs to use.
	&wfs
	nwn = 15, nwp = 15
	/
	# Inputfilename  = file from which to continue the calculation
	# Outputfilename = .wf file to write after the end of the calculation. 
	# init signals the code to perform its own initialization.
	&IO
	InputFilename='init'
	Outputfilename='tant.wf'
	/
	&MomentParam
	/
	&Cranking
	/
	EOF

For more information, please take look in the `tantalus_examples/` folder, which contains its own readme with more deails on the calculations.

A manual detailing (some of) the more advanced runtime options of Tantalus will be available soon<sup>TM</sup>.

--------
### Configuration files

Configuration files in the `configs/` folder dictate the mean-field code produced at the end of the day by Hephaestos. The default example is


	# Configuration file for Tantalus compilation
	# Functional : NLO Skyrme-type (think SLy4/5/6) with optional Jmn terms
	FUNC_FILE = 'NLO.func'
	# Symmetries : EV8-style
	SYMSTRING = 'Rz,T,P,STy'
	REDUCE    = [1,1,1]
	# Read symmetries: EV8-style
	INSYM     = 'Rz,T,P,STy'
	INREDUCE  = [1,1,1]

* **FUNC_FILE** : indicates the functional file that will be implemented, from the `functionals/` folder.
* **SYM_STRING**: the generators of the single-particle symmetry group for the calculation, i.e. the symmetries which give us *single-particle* symmetry relations (quantum numbers in the case of linear operators). 
* **REDUCE**: A python list of 0 or 1's, indicating which Cartesian axis we want to reduce the computation of. I.e. 1 means only half of the axis is stored, 0 means the axis is used in its entirety.
* **IN_SYM** : symmetry string in the same format  as SYM_STRING. Dictates the single-particle symmetries that will be expected on a .wf file read on input. 
* **IN_REDUCE**: reduction string of the same format as REDUCE. Dictates the reduction of Cartesian axes expected on a .wf file read on input. 

Note that Tantalus is **always prepared to read .wf files from a previous calculation with the same .exe.** The IN_SYM and IN_REDUCE variables only affect what other type of .wfs file can be read, i.e. what to do if a .wf is read from a calculation with **DIFFERENT** symmetry choices.