################################################################################
# Example runscript for Tantalus. 
#
# Here, we demonstrate a cranking calculation.
#
# Run: 
#     (1) Calculate a deformed configuration for 36Ar
#         Note: prior experience has shown this particular initialisation 
#               of the constraints ends up in an oblate configuration 
#               with y as the short axis. 
#   
################################################################################

exe='Tantalus.default.exe'
execdir='../exec'
paramloc='../parameterizations/'
param="SLy5s1"
iterations=300

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
echo " --> Deformed starting configuration"
outfile="Tant.Ar36.out"
#-------------------------------------------------------------------------------
cat << EOF > tant.Ar36.data
&nucleus
neutrons=18, protons=18
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type="HF"
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
Outputfilename='tant.Ar36.wf'
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=25.0
iteration=10
/
&Cranking
/
EOF
#-------------------------------------------------------------------------------
# Running the code
./$exe < tant.Ar36.data > $outfile
mv $outfile ../out/STDOUT/

################################################################################
exe='Tantalus.NLO-T.exe'
cd ..
cp $execdir/$exe work/
cd  work/
#-------------------------------------------------------------------------------
mv tant.Ar36.wf tant.in
for omz in 0.1 0.2 0.3 0.4 0.5
do

echo " --> Running a cranking calculation, om_z = $omz"
outfile="Tant.Ar36.om_z=$omz.out"
#-------------------------------------------------------------------------------

cat << EOF > tant.Ar36.data
&nucleus
neutrons=18, protons=18
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$param"
/
&pairing
Type="HF"
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
InputFilename='tant.in'
Outputfilename="tant.Ar36.om=$omz.wf"
allowtransform=.true.
/
&MomentParam
/
&MomentConstraint
/
&Cranking
omegaz = $omz
/
EOF

./$exe < tant.Ar36.data > $outfile
mv tant.Ar36.om=$omz.wf tant.in

done
mv *.out ../out/STDOUT
rm *.exe *.param *.dat

