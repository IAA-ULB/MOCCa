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
#  This particular runscript illustrates taking the input from a previous 
#  calculation to initialize a second one.
#
################################################################################

exe='Tantalus.default.exe'
execdir='../exec'
paramloc='../parameterizations/'
param="SLy5s1"

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

outfile="Tant.Ca40.out"
echo " --> Starting a calculation for 40Ca"

#-------------------------------------------------------------------------------
cat << EOF > tant.Ca40.data
&nucleus
neutrons=20, protons=20
/
&mesh
nx=14, ny=14, nz=14, dx=0.8
/
&func
name_param='SLy5s1'
/
&pairing
Type="HFB"
/
&evolution
maxiter=250
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
/
&IO
InputFilename='init'
Outputfilename='tant.Ca40.wf'
BXLFIT='continuing.'
/
&MomentParam
/
&Cranking
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Ca40.data > $outfile

#Cleaning up
mv $outfile     ../out/STDOUT
mv continuing.* ../out/summary
################################################################################

outfile="Tant.Ca42.out"
echo "  --> Starting a calculation for 42Ca from the previous one"

#-------------------------------------------------------------------------------
cat << EOF > tant.Ca42.data
&nucleus
neutrons=22, protons=20
/
&mesh
nx=14, ny=14, nz=14, dx=0.8
/
&func
name_param='SLy5s1'
/
&pairing
Type="BCS"
/
&evolution
maxiter=250
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
/
&IO
InputFilename='tant.Ca40.wf'
Outputfilename='tant.Ca42.wf'
BXLFIT='continuing.'
/
&MomentParam
/
&Cranking
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Ca42.data > $outfile

#Cleaning up
rm tant.*.wf  
mv $outfile     ../out/STDOUT
mv continuing.* ../out/summary
#-------------------------------------------------------------------------------
rm *.exe *.data *.param
