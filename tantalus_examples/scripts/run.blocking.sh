################################################################################
# Example runscript for Tantalus. 
#
# Here, we demonstrate the possibility for various options for qp blocking
#
# Run: 
#     (1) Converges a false vacuum for Mg25
#     (2) Uses the previous run to initialize an EFA calculation for Mg25
#     (3) Uses that final run to perform a time-reversal broken calculation, 
#         still blocking EFA-style
#     (4) Finally, transfer the EFA-blocking to a real blocked calculation
#
################################################################################

exe='Tantalus.BXL.exe'
execdir='../exec'
paramloc='../parameterizations/'
param="BSkG1"
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
outfile="Tant.Mg25.FV.out"
#-------------------------------------------------------------------------------
cat << EOF > tant.Mg25.data
&nucleus
neutrons=13, protons=12
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
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
nwn = 20, nwp = 20
/
&IO
InputFilename='init'
Outputfilename='tant.Mg25.FV.wf'
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=8.0
iteration=10
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Mg25.data > $outfile
mv $outfile ../out/STDOUT/

#-------------------------------------------------------------------------------
echo " --> EFA calculation"
outfile="Tant.Mg25.EFA.out"

cat << EOF > tant.Mg25.data
&nucleus
neutrons=13, protons=12
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type="HFB"
blocktype=4
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
nwn = 20, nwp = 20
/
&IO
InputFilename='tant.Mg25.FV.wf'
Outputfilename='tant.Mg25.EFA.wf'
/
&MomentParam
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Mg25.data > $outfile
mv $outfile ../out/STDOUT/
#-------------------------------------------------------------------------------
echo " --> EFA blocking without T calculation"
outfile="Tant.Mg25.EFA.T.out"

exe='Tantalus.BXL-T.exe'
execdir='../exec'
cd ..
cp $execdir/$exe             work/
cd  work

cat << EOF > tant.Mg25.data
&nucleus
neutrons=13, protons=12
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type="HFB"
blocktype=4
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
InputFilename='tant.Mg25.EFA.wf'
Outputfilename='tant.Mg25.EFA.T.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Mg25.data > $outfile
mv $outfile ../out/STDOUT/
#-------------------------------------------------------------------------------
echo " --> Full blocking"
outfile="Tant.Mg25.block.out"

cat << EOF > tant.Mg25.data
&nucleus
neutrons=13, protons=12
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
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
InputFilename='tant.Mg25.EFA.T.wf'
Outputfilename='tant.Mg25.block.wf'
/
&MomentParam
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Mg25.data > $outfile
mv $outfile ../out/STDOUT/
rm *.wf *.exe *.param
