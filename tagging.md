# Blocking with tagging

## Generating a tagging spwf
To generate a single-particle wavefunction to use as a tag, you should use the auxiliary gen_nilsson program. This can be compiled by typing 

> make gen_nilsson 

in the main directory. Standard options regarding compilers apply. 

You can run the resulting executable (exec/gen_nilsson.exe) with the following namelist input

> nx = 12         !
> ny = 12         ! mesh parameters 
> nz = 12         !
> dx = 1.0        !
> neutrons = 58   ! particle numbers for the Nilsson Hamiltonian
> protons  = 41   ! 
> nwn = 150       ! total number of neutron spwfs
> nwp = 150       ! total number of proton spwfs
> osc_x = 0.2     ! oscillator frequencies of the Nilsson Hamiltonian
> osc_y = 0.2     !
> osc_z = 0.15    !
> fname='spwf.dat'! Filename to write to
> selection=165   ! Selected index, i.e. the index of the spwf that will be written to file.

When this data is used, you will end up with 'spwf.dat' file (or whatever you specified as fname) that 
can be used in a future MOCCa run as a tag.

Note: Gen_nilsson relies on the same routines as the MOCCa initialisation routines; i.e. if you specify the same input, the spwfs at iteration 0 of a MOCCa run are IDENTICAL to those you could get here?  

## Tagging with MOCCa

Assuming you have constructed an appropriate tagging file, you can tell MOCCa to use it by specifying the following options:

> &pairing
> Type          = HFB
> tag_spwf_file = 'spwf.dat'  ! filename of the tagging file
> blocknumber   =1            ! number of qp's to block
> blocktype=7                 ! signal that tagging is the blocking strategy
> /
> &indices
> /

Note that the &indices/ namelist can remain empty, in contrast to all other blocking options!

## Limitations of the implementation

The tagging strategy currently only works for 
- HFB calculations. HF and HF+BCS is not supported.
- Time-reversal broken calculations. Time-reversal conserved calculations relying on the equal
  filling approximations are not supported.
