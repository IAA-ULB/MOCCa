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
#  This particular runscript illustrates a calculation for Ca48 with SN2LO1.
#
################################################################################

exe='Tantalus.N2LO.func.exe'
execdir='../exec'
param="SN2LO1"
paramloc="../parameterizations/"
outfile='Tant.N2LO.out'

#Create storage directories
if [ ! -d "out/" ]; then
    mkdir out
fi
if [ ! -d "out/STDOUT" ]; then
   mkdir out/STDOUT
fi
if [ ! -d "out/summary" ]; then
   mkdir out/summary
fi
if [ ! -d "work/" ]; then
  mkdir work
fi

cp $execdir/$exe             work/
cp $paramloc/"$param.param"  work/

cd work

#-------------------------------------------------------------------------------
# Creating the runtime data
cat << EOF > tant.data
&nucleus
neutrons=28, protons=20
/
# Parameters of the Lagrange mesh. 
&mesh
nx=14, ny=14, nz=14, dx=0.8
/
# The code will look, on a file forces.param, for the parameterization with 
# this name.
&func
name_param="$param"
/
# Options for the pairing.
&pairing
/
# maxiter = Maximum number of iterations to be performed
&evolution
maxiter=250
/
&scfiteration
/
# Number of neutron (nwn) and proton (nwp) spwfs to use.
&wfs
nwn = 35, nwp = 35
/
# Inputfilename  = file from which to continue the calculation
# Outputfilename = .wf file to write after the end of the calculation. 
# init signals the code to perform its own initialization.
&IO
InputFilename='init'
Outputfilename='tant.wf'
BXLFIT='N2LO.'
/
&MomentParam
/
&Cranking
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.data > $outfile

#Cleaning up
rm tant.wf  
mv $outfile  ../out/STDOUT
mv N2LO.*    ../out/summary
rm *.exe *.data *.param
#-------------------------------------------------------------------------------
