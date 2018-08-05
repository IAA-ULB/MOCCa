

# Tantalus & Hephaestos
> Copyright W. Ryssens, P.H. Heenen & M. Bender
> Sunday, 05. August 2018 09:02PM 

----------
### An introduction to Greek mythology

Or, perhaps more appropriate
> I have been given access to this repository, and have no clue where to start?!

__Tantalus__ is a mean-field code, representing the single-particle wavefunctions on a Lagrange mesh, written in FORTRAN.
__Hephaestos__ is a (collection of) Python modules/scripts that partially writes (parts of) the FORTRAN source code for Tantalus. 


Hephaestos, is in fact a code generator or preprocessor for Tantalus. For different choices regarding types of functionals and symmetries, Hephaestos generates a version of Tantalus that is able to perform calculations for you, subject to the input you provided Hephaestos.

-------
### File structure

Files in the repository's main folder are

*   `Hephaestos.py`: Main python script of Hephaestos.
*  `Makefile`: Makefile to compile Tantalus. 

Source code folders are

*  `src_heph/` :  contains modules for Hephaestos (in python). 
*  `src_orig/` :  contains source code files (in FORTRAN) of Tantalus, to be completed by Hephaestos.
*  `src/` :   directory of the compileable Tantalus source code, where Hephaestos places the processed .f90 files.

Information on the various functionals are

*  `functionals/` : Contains the definition of various forms of functionals (`.func` files), as input for Hephaestos.
*  `parameterizations/` : Contains various parameterization file (`.param` files), each linked to a particular `.func` file.  

Auxiliary folders are

* `exec/` :  The Makefile places the compiled executables here.
* `mod/` : Auxiliary folder where the Makefile places the compiled `.mod` files.
* `obj/` :  Auxiliary folder fo the compiled `.obj` files. 

Other

* `tantalus_examples/` : contains example-scripts for running Tantalus.

-------
### How do I compile the code?

The easiest way is to type

`make`

By default, this 

1) Runs Hephaestos, for a defined type of functional.
2) The processed source code files are placed in `̀src/`.
2) Compiles the processed source code
3) The compiled exe is placed exe into `exec/`.

compiles Tantalus for use with an NLO functional in its maximally symmetric run-mode. 

Currently defined options to pass into the Makefile are

* `CXX`: your preferred choice of compiler. Default = `gfortran`.
* `CXXFLAGS`:  options for your preferred compiler.  
* `FUNC`:  name of the functional file to take as input for Hephaestos. Default = `NLO.func` 

__Things do check if the code does not want to compile__

* Are the `src/`,`mod/` and `obj/` directories present in your directory tree?

-------
### How do I run Tantalus? 

The minimal STDIN input for a calculation is given

For more information, please take look in the `tantalus_examples/` folder. 

A manual detailing (some of) the more advanced runtime options of Tantalus will be available soon<sup>TM</sup>.

-------
### How do I run Hephaestos?

Answer:	
		
		python Hephaestos.py NLO.func

where `NLO.func` is the name of the functional file in the `functionals` folder, and can be replaced by the functional definition file of your choice.  (Specifying a functional file is __not optional__!) The code will process files from the `src_orig/ ` directory and place processed files in the `src` directory. 

Separate running of Hephaestos is not immediately useful, as after that Tantalus still needs to be compiled. 