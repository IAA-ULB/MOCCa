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
# This is an example script showing the calculation of the magnetic dipole 
# moment of the 7/2- ground-state in Sc41 with the SLy5s1 EDF.
#
################################################################################

exe='Tantalus.default.exe'
exe_T='Tantalus.NLO-T.exe'
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
cp $execdir/$exe     work/
cp $execdir/$exe_T   work/

param=SLy5s1
cp ../parameterizations/$param.param  work/$param.param

cd work

outfile="Tant.Sc41.$param.FV.out"

#-------------------------------------------------------------------------------
# Creating the runtime data for a false vacuum
cat << EOF > tant.data
&nucleus
neutrons=20, protons=21
/
&mesh
nx=12, ny=12, nz=12, dx=0.8
/
&func
name_param="$param"
/
# Options for the pairing.
&pairing
type='HFB'
/
&evolution
maxiter=100
printiter=10
/
&scfiteration
/
&wfs
nwn = 30, nwp = 30
/
&IO
InputFilename='init'
Outputfilename='tant.FV.wf'
/
&MomentParam
/
&MomentConstraint
/
&Cranking
/
EOF

#-------------------------------------------------------------------------------
# Running the code for a false vacuum
echo "Running a false vacuum calculation first."
./$exe < tant.data > $outfile
mv $outfile ../out/STDOUT

#-------------------------------------------------------------------------------
# Creating the runtime data for a fully blocked calculation
cat << EOF > tant.data
&nucleus
neutrons=20, protons=21
/
&mesh
nx=12, ny=12, nz=12, dx=0.8
/
&func
name_param="$param"
/
# Options for the pairing.
&pairing
type='HFB'
blocknumber=1
blocktype=1    
/
&indices
blockindices=115
/
&evolution
maxiter=100
printiter=10
/
&scfiteration
/
# Double the # of spwfs for use in a time-reversal broken calculation.
&wfs
nwn = 60, nwp = 60
/
# AllowTransform=.true., such that the code can break the T-symmetry
&IO
InputFilename='tant.FV.wf'
Outputfilename='tant.blocked.wf'
allowtransform=.true.
/
&MomentParam
/
&MomentConstraint
/
&Cranking
/
EOF

echo "Running the blocked calculation."
outfile="Tant.Sc41.$param.blocked.out"
./$exe_T < tant.data > $outfile


#Cleaning up
rm tant.*.wf 
rm *.exe *.data 
#-------------------------------------------------------------------------------
