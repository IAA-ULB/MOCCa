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
# This particular script tests a calculation of 208Pb for older Brussels 
# interactions with and without t4/t5 terms. 
#
#  BSk27
#  BSk29
#  BSk31
################################################################################

exe='Tantalus.BXL.exe'
execdir='../exec'
paramloc='../parameterizations/'

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

cd work

for param in BSk27 BSk29 BSk31
do

cp ../$paramloc/"$param.param"  .

echo " --> Calculation for Pb208 with $param"
outfile="Tant.Pb208.$param.out"

#-------------------------------------------------------------------------------
cat << EOF > tant.data
&nucleus
neutrons=126, protons=82
pairing_prec=1e-10
/
&mesh
nx=18, ny=18, nz=18, dx=0.7
/
# The SLy4 parameterization has pairing strength and type defined in its 
# .param file.
&func
name_param='$param'
/
&pairing
Type="HF"
/
&evolution
maxiter=200
/
&scfiteration
/
&wfs
nwn = 70, nwp = 50
osc_freq = 0.2,0.2,0.2
/
&IO
InputFilename='init'
Outputfilename='tant.wf'
BXLFIT='Pb208.'
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
mv $outfile ../out/STDOUT
mv Pb208.*  ../out/summary
done
#-------------------------------------------------------------------------------
rm *.exe
rm *.data

