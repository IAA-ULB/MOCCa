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
#  This particular runscript illustrates a two calculations, employing BCS and
#  HFB ansatzes for 24Mg in a limited box.
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

for type in BCS HFB
do

echo "Calculation with pairingtype=$type"
outfile="Tant.$type.out"


#-------------------------------------------------------------------------------
cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
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
maxiter=100
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
/
&IO
InputFilename='init'
Outputfilename='tant.wf'
/
&MomentParam
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.data > $outfile

#Cleaning up
mv tant.wf  ../wf
mv $outfile ../out
#-------------------------------------------------------------------------------
done
rm *.exe
rm *.data
