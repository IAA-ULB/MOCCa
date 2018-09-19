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
#  This particular runscript illustrates a constrained calculation, obtaining
#  a triaxial configuration for 20Ne.
#
################################################################################

exe='Tantalus.NLO.func.exe'
execdir='../exec'
param='../parameterizations/SLy5s1.param'

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
cp $param          work/forces.param

cd work


outfile="Tant.Ca40.out"
echo "Starting a calculation for 40Ca"

#-------------------------------------------------------------------------------
cat << EOF > tant.Ca40.data
&nucleus
neutrons=20, protons=20
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy5s1'
/
&pairing
Type="BCS"
/
&evolution
maxiter=100
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
/
&IO
InputFilename='init'
Outputfilename='tant.Ca40.wf'
/
&MomentParam
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Ca40.data > $outfile

#Cleaning up
mv $outfile ../out

################################################################################

outfile="Tant.Ca42.out"
echo "Starting a calculation for 42Ca from the previous one"

#-------------------------------------------------------------------------------
cat << EOF > tant.Ca42.data
&nucleus
neutrons=22, protons=20
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy5s1'
/
&pairing
Type="BCS"
/
&evolution
maxiter=100
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
/
&IO
InputFilename='tant.Ca40.wf'
Outputfilename='tant.Ca42.wf'
/
&MomentParam
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Ca42.data > $outfile

#Cleaning up
mv tant.*.wf  ../wf
mv $outfile ../out

#-------------------------------------------------------------------------------
rm *.exe
#rm *.data
rm forces.param
