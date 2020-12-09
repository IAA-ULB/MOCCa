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
#  This particular runscript illustrates two calculations, employing BCS and
#  HFB ansatzes for Ca44 in a limited box.
#
################################################################################

exe='Tantalus.default.exe'
execdir='../exec'
param="SLy5s1"
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
cp $paramloc/"$param.param"  work/

cd work

for type in BCS HFB
do

echo "Calculation with pairingtype=$type"
outfile="Tant.$type.out"


#-------------------------------------------------------------------------------
cat << EOF > tant.data
&nucleus
neutrons=24, protons=20
/
&mesh
nx=14, ny=14, nz=14, dx=0.8
/
# The SLy5s1 has pairing strength and type defined in its .param file.
&func
name_param='SLy5s1'
/
# Start with BCS pairing.
&pairing
Type="$type"
/
&evolution
maxiter=300
/
&scfiteration
/
&wfs
nwn = 30, nwp = 30
/
&IO
InputFilename='init'
Outputfilename='tant.wf'
BXLFIT="pairing.$type."
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
mv $outfile  ../out/STDOUT
mv pairing.* ../out/summary
#-------------------------------------------------------------------------------
done
rm *.exe *.data *.param

