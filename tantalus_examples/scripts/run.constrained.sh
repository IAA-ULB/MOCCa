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
param="SLy5s1"
paramloc='../parameterizations'

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

cp $execdir/$exe            work/
cp $paramloc/"$param.param" work/

cd work

for Q20 in `seq 10 10 30`
do

outfile="Tant.Q20=$Q20.out"
echo "Calculating configuration with Q20=$Q20 fm^2"

#-------------------------------------------------------------------------------
cat << EOF > tant.data
&nucleus
neutrons=10, protons=10
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
maxiter=300
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
/
&IO
InputFilename='init'
Outputfilename='tant.wf'
BXLFIT="constrained.Q20=$Q20."
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=0
Constraint=$Q20
/
EOF

#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.data > $outfile

#Cleaning up
rm tant.wf  
mv $outfile      ../out/STDOUT
mv constrained.* ../out/summary
#-------------------------------------------------------------------------------
done
rm *.exe *.data *.param

