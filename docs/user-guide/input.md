# Running MOCCa
A boilerplate MOCCa calculation requires that you provide at least:

* [a **.in file**: specifying all runtime parameters](structure.md).
* [a **.param file**: the details of the EDF parameterisation](param.md).

You can optionally provide the following files to warmstart the code:

* [a **.wf file**: a wavefunction file obtained in a previous MOCCa calculation ](wffile.md).
* [a **.pot file**:  mean-field potentials obtained from a simpler type of calculation ](potfile.md).


Actually executing the code can then be as simple as typing

    ./MOCCa.exe data.in 
    
where you have to make sure that specification of the other input files in 'data.in' reference locations that the executable can access.

Note that it is also possible to pipe in the data in the .in file as follows:

    ./MOCCa.exe < data.in
