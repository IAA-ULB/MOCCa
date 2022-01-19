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

exe='Tantalus.BXL-T.exe'
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
cp ../parameterizations/BSkG1.param         work/BSkG1.param
cp ../parameterizations/BSkG2.param         work/BSkG2.param

cd work

for param in BSkG1 BSkG2
do

echo " --> Running $param"
outfile="Tant.$param.T.out"

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
type='HFB'
/
# maxiter = Maximum number of iterations to be performed
&evolution
maxiter=100
printiter=10
/
&scfiteration
/
# Number of neutron (nwn) and proton (nwp) spwfs to use.
&wfs
nwn = 60, nwp = 60
/
# Inputfilename  = file from which to continue the calculation
# Outputfilename = .wf file to write after the end of the calculation. 
# init signals the code to perform its own initialization.
&IO
InputFilename='init'
Outputfilename='tant.wf'
BXLFIT="$param.T."
allowtransform=.true.
/
&MomentParam
MoreConstraints=.false.
/
&MomentConstraint
l=2
m=0
constraint=10.0
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=0
/
&Cranking
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.data > $outfile
mv $outfile ../out/STDOUT
mv $param.* ../out/summary

done

#Cleaning up
rm tant.wf 
rm *.exe *.data 
#-------------------------------------------------------------------------------
