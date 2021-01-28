################################################################################
# Example runscript for Tantalus. 
#
# Here, we demonstrate the possibility for more advanced forms of wavefunction
# input. 
#
# Run: 
#     (1) Create a .wf file from scratch, including time-reversal.
#     (2) Start from that .wf file, and add a few extra spwfs. 
#     (3) Using the extended .wf file, break time-reversal
#     (4) Extend that .wf file again with a few extra spwfs
#
################################################################################

exe='Tantalus.default.exe'
execdir='../exec'
paramloc='../parameterizations/'
param="SLy5s1"
iterations=100

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
echo " --> Creating a .wf file from scratch"
outfile="Tant.O16.out"
#-------------------------------------------------------------------------------
cat << EOF > tant.O16.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy5s1'
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
nwn = 10, nwp = 10
/
&IO
InputFilename='init'
Outputfilename='tant.O16.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.O16.data > $outfile

################################################################################
echo " --> Adding some extra spwfs on input"
outfile="Tant.O16.extra.out"
#-------------------------------------------------------------------------------
cat << EOF > tant.O16.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy5s1'
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
nwn = 14, nwp = 14
/
&IO
InputFilename='tant.O16.wf'
Outputfilename='tant.O16.extra.wf'
AllowTransform=.true.
extraspwfs= 2,0,2,0,2,0,2,0
/
&MomentParam
/
&Cranking
/
EOF

# Running the code
./$exe < tant.O16.data > $outfile
################################################################################
echo " --> Breaking Time-reversal"
cd ..
exe='Tantalus.NLO-T.exe'
cp $execdir/$exe work/
cd  work

outfile="Tant.O16.T.out"
#-------------------------------------------------------------------------------
cat << EOF > tant.O16.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy5s1'
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
nwn = 28, nwp = 28
/
&IO
InputFilename='tant.O16.extra.wf'
Outputfilename='tant.O16.T.wf'
AllowTransform=.true.
/
&MomentParam
/
&Cranking
/
EOF

# Running the code
./$exe < tant.O16.data > $outfile
################################################################################
echo " --> Adding spwfs to the time-reversal-broken case"

outfile="Tant.O16.T.extra.out"
#-------------------------------------------------------------------------------
cat << EOF > tant.O16.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy5s1'
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
nwn = 36, nwp = 36
/
&IO
InputFilename='tant.O16.T.wf'
Outputfilename='tant.O16.T.extra.wf'
AllowTransform=.true.
extraspwfs = 2,2,2,2,2,2,2,2
/
&MomentParam
/
&Cranking
/
EOF

# Running the code
./$exe < tant.O16.data > $outfile
mv *.out ../out/STDOUT/
rm *.exe *.data *.param
rm *.wf
