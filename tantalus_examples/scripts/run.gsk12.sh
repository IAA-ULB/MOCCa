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
# This is an example script using the GSk1 and GSk2 parameterizations for Ca48.
#
################################################################################

exe='Tantalus.DD-switch.exe'
execdir='../exec'

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
cp $execdir/$exe   work/
cp ../parameterizations/gsk1.param         work/gsk1.param
cp ../parameterizations/gsk2.param         work/gsk2.param

cd work

for param in gsk1 gsk2
do

echo " --> Running $param"
outfile="Tant.$param.out"

#-------------------------------------------------------------------------------
# Creating the runtime data
cat << EOF > tant.data
&nucleus
neutrons=28, protons=20
/
# Parameters of the Lagrange mesh. 
&mesh
nx=16, ny=16, nz=16, dx=0.8
/
# The code will look, on a file forces.param, for the parameterization with 
# this name.
&func
name_param="$param"
/
# Options for the pairing.
&pairing
type='HFB'
/
# maxiter = Maximum number of iterations to be performed
&evolution
maxiter=100
/
&scfiteration
/
# Number of neutron (nwn) and proton (nwp) spwfs to use.
&wfs
nwn = 30, nwp = 30
/
# Inputfilename  = file from which to continue the calculation
# Outputfilename = .wf file to write after the end of the calculation. 
# init signals the code to perform its own initialization.
&IO
InputFilename='init'
Outputfilename='tant.wf'
BXLFIT="$param."
/
&MomentParam
/
&Cranking
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.data > $outfile
mv $outfile ../out/STDOUT
mv $param.z* ../out/summary

done

#Cleaning up
rm tant.wf 
rm *.exe *.data *.param
#-------------------------------------------------------------------------------
