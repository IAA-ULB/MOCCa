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
# This particular runscript tests a calculation of 208Pb in realistic conditions.
#
################################################################################

exe='Tantalus.NLO.func.exe'
execdir='../exec'
paramloc='../parameterizations/'
param='SLy4'

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

echo "Realistic calculation for Pb208"
outfile="Tant.Pb208.out"

#-------------------------------------------------------------------------------
cat << EOF > tant.data
&nucleus
neutrons=126, protons=82
pairing_prec=1e-10
/
&mesh
nx=16, ny=16, nz=16, dx=0.8
/
# The SLy4 parameterization has pairing strength and type defined in its 
# .param file.
&func
name_param='SLy4'
/
&pairing
Type="HFB"
/
&evolution
maxiter=200
/
&scfiteration
/
&wfs
nwn = 70, nwp = 50
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
#-------------------------------------------------------------------------------
rm *.exe
rm *.data
