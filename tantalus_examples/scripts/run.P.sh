################################################################################
#
# This is an example script illustrating the use of the code for parity-broken
# calculations. 
#
################################################################################

exe='Tantalus.DD-switch.exe'
exeb='Tantalus.DD-switch-P.exe'
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
cp $execdir/$exeb   work/

cp ../parameterizations/gsk2.param         work/gsk2.param

cd work

param=gsk2
outfile="Tant.$param.symmetric.out"

#-------------------------------------------------------------------------------
# Creating an initial (symmetric) wavefunction file
cat << EOF > tant.data
&nucleus
neutrons=10, protons=10
/
# Parameters of the Lagrange mesh. 
&mesh
nx=12, ny=12, nz=12, dx=0.8
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
maxiter=200
printiter=10
/
&scfiteration
/
# Number of neutron (nwn) and proton (nwp) spwfs to use.
&wfs
nwn = 20, nwp = 20
/
# Inputfilename  = file from which to continue the calculation
# Outputfilename = .wf file to write after the end of the calculation. 
# init signals the code to perform its own initialization.
&IO
InputFilename='init'
Outputfilename='tant.wf'
allowtransform=.true.
/
&MomentParam
MoreConstraints=.false.
/
&Cranking
/
EOF

#-------------------------------------------------------------------------------
# Running the code
echo "-> Running a symmetric calculation"
./$exe < tant.data > $outfile
mv $outfile ../out/STDOUT
#-------------------------------------------------------------------------------

outfile="Tant.$param.Q30=30.out"

cat << EOF > tant.data
&nucleus
neutrons=10, protons=10
/
# Parameters of the Lagrange mesh. 
&mesh
nx=12, ny=12, nz=24, dx=0.8
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
maxiter=200
printiter=10
/
&scfiteration
/
# Number of neutron (nwn) and proton (nwp) spwfs to use.
&wfs
nwn = 20, nwp = 20
/
# Inputfilename  = file from which to continue the calculation
# Outputfilename = .wf file to write after the end of the calculation. 
# init signals the code to perform its own initialization.
&IO
InputFilename='tant.wf'
Outputfilename='tant.wf'
allowtransform=.true.
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=3
m=0
constraint=50.0
/
&Cranking
/
EOF

echo "-> Breaking parity"
./$exeb < tant.data > $outfile
mv $outfile ../out/STDOUT

#-------------------------------------------------------------------------------
#Cleaning up
rm tant.wf 
rm *.exe *.data  *.param
#-------------------------------------------------------------------------------
