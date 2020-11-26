################################################################################
# Example runscript for Tantalus. 
#
# For succesful execution you need:
#
#   a) An exe, compiled with a given .func file. 
#   b) A .param file, that corresponds to the .func file used for compilation
#   c) A data file, passed into Tantalus via STDIN. 
#   d) A tantalus .wf file from a previous calculation, if you do not want to 
#      start from scratch. 
#-------------------------------------------------------------------------------
#
#  This particular runscript illustrates the most basic calculation.
#  Initializing from a set of Nilsson orbitals (inputfilename='init'), the 
#  code performs a Hartree-Fock optimization for 16O, with the SLy4 
#  parameterization of the NLO functional, in a limited box of (12,12,12,1.0)
#  points with all symmetries possible conserved.
#
################################################################################

exe='Tantalus.default.exe'
execdir='../exec'
param='../parameterizations/SLy4.param'
outfile='Tant.minimal.out'

#Create storage directories
if [ ! -d "out/" ]; then
  mkdir out
fi
if [ ! -d "wf/" ]; then
  mkdir wf
fi
if [ ! -d "work/" ]; then
  mkdir work
fi

cp $execdir/$exe   work/
cp $param          work/SLy4.param

cd work

#-------------------------------------------------------------------------------
# Creating the runtime data
cat << EOF > tant.data
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

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.data | tee $outfile

#Cleaning up
mv tant.wf  ../wf
mv $outfile ../out
rm *.exe
rm *.data
#-------------------------------------------------------------------------------
