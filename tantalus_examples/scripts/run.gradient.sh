################################################################################
# Example runscript for Tantalus. 
#
# Here, we demonstrate the use of the gradient solver for the HFB problem.
# 1) Converge a false vacuum for O19
# 2) Block the lowest n+ state and try to converge the direct case
# 3) Block the lowest n+ state and converge with the gradient solver
################################################################################

exe='Tantalus.DD-switch-T.exe'
execdir='../exec'
paramloc='../parameterizations/'
param="gsk2"
iterations=200

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

################################################################################
echo " --> False vacuum calculation"
outfile="Tant.O19.FV.out"
#-------------------------------------------------------------------------------
cat << EOF > tant.O19.data
&nucleus
neutrons=9, protons=8
/
&mesh
nx=10, ny=10, nz=10, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type="HFB"
/
&evolution
maxiter=$iterations 
/
&scfiteration
/
&wfs
nwn = 40, nwp = 40
/
&IO
InputFilename='init'
Outputfilename='tant.O19.FV.wf'
allowtransform=.true.
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=8.0
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=0.0
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
## Running the code
./$exe < tant.O19.data > $outfile
#mv $outfile ../out/STDOUT/
#-------------------------------------------------------------------------------
echo " --> Direct blocking"
outfile="Tant.O19.direct.out"

cat << EOF > tant.O19.data
&nucleus
neutrons=9, protons=8
/
&mesh
nx=10, ny=10, nz=10, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type="HFB"
blocktype=2
blocknumber=1
/
&indices
blocklowest='n+'
/
&evolution
maxiter=$iterations 
/
&scfiteration
/
&wfs
nwn = 40, nwp = 40
/
&IO
InputFilename='tant.O19.FV.wf'
Outputfilename='tant.O19.direct.wf'
allowtransform=.true.
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=8.0
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=0.0
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.O19.data > $outfile
mv $outfile ../out/STDOUT/

#-------------------------------------------------------------------------------
echo " --> Gradient blocking"
outfile="Tant.O19.grad.out"

cat << EOF > tant.O19.data
&nucleus
neutrons=9, protons=8
/
&mesh
nx=10, ny=10, nz=10, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type="HFB"
blocktype=2
blocknumber=1
pairingscheme=1
bogofromfile=.false.
/
&indices
blocklowest='n+'
/
&evolution
maxiter=$iterations 
/
&scfiteration
/
&wfs
nwn = 40, nwp = 40
/
&IO
InputFilename='tant.O19.FV.wf'
Outputfilename='tant.O19.gradient.wf'
allowtransform=.true.
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=8.0
MoreConstraints=.true.
multfromfile=.true.
/
&MomentConstraint
l=2
m=2
constraint=0.0
multfromfile=.true.
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.O19.data > $outfile
mv $outfile ../out/STDOUT/
rm *.wf *.exe *.param
