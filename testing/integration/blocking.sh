#-------------------------------------------------------------------------------
# Test the blocking functionality for 55Cr
#
# This script tests the 4-step blocking workflow:
#  1. False vacuum calculation
#  2. EFA calculation - using false vacuum as input
#  3. EFA blocking without T calculation - using previous calculation as input
#  4. Full blocking
#
# and compares them to known results.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash blocking.sh exec exec_T
#
# where
# - exec specifies the suffix of the executable to be used when T is conserved
#            Example: "BXL" for "MOCCa.BXL.exe".
#
# - exec_T specifies the suffix of the executable to be used when T is broken
#            Example: "BXL-T" for "MOCCa.BXL-T.exe".
#
# Dependencies: BSkG1 parameterization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit     : e579b30ca1604b89dd2ccec045f73b811dca49b7
#-------------------------------------------------------------------------------

# KNOWN RESULTS TO COMPARE AGAINST
# Add your expected values here for each of the 4 calculations

# Step 1: False vacuum calculation
refE_FV=-479.130 # Total energy for false vacuum in MeV
refB20_FV=0.000  # Quadrupole deformation beta_20
refB22_FV=0.000  # Quadrupole deformation beta_22

# Step 2: EFA calculation
refE_EFA=-479.205 # Total energy for EFA in MeV
refB20_EFA=0.090  # Quadrupole deformation beta_20
refB22_EFA=0.000  # Quadrupole deformation beta_22

# Step 3: EFA blocking without T calculation
refE_EFA_T=-479.208 # Total energy for EFA-T in MeV
refB20_EFA_T=0.090  # Quadrupole deformation beta_20
refB22_EFA_T=0.000  # Quadrupole deformation beta_22

# Step 4: Full blocking
refE_block=-478.890 # Total energy for full blocking in MeV
refB20_block=0.083  # Quadrupole deformation beta_20
refB22_block=0.000  # Quadrupole deformation beta_22

# Tolerances for comparison
E_tol=0.005 # Energy
B_tol=0.02  # Deformation

#set -v
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

param="BSkG1"
iterations=1000

#-------------------------------------------------------------------------------
# Set up a single work directory for all 4 steps
# We need to manually set up since we need to share wavefunctions between steps
exesuffix=$1
exesuffix_T=$2
exebase="MOCCa.$exesuffix.exe"

if [ ! -d "work/" ]; then
  mkdir work
fi
if [ ! -d "logs/" ]; then
  mkdir logs
fi
if [ ! -d "../logs/" ]; then
  mkdir ../logs
fi

cp ../../exec/$exebase work/
cp ../../exec/MOCCa.${exesuffix_T}.exe work/ 2>/dev/null || true
cp ../../parameterizations/$param.param work/

cd work

# Clean up any previous wavefunctions
rm -f MOCCa.Cr55.*.wf

echo " --> Step 1: False vacuum calculation"

cat <<EOF >MOCCa.Cr55.data
&nucleus
neutrons=31, protons=24
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
nwn = 30, nwp = 30
/
&IO
InputFilename='init'
Outputfilename='MOCCa.Cr55.FV.wf'
/
&MomentParam
MoreConstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=80.0
iteration=10
/
&Cranking
/
EOF

# Running the code
./$exebase <MOCCa.Cr55.data >../logs/blocking_FV.out
mocca_check_FV=$?
outfile="../logs/blocking_FV.out"

# Check results for Step 1
E_FV=$(get_total_energy_stdout $outfile)
B20_FV=$(get_B20_stdout $outfile)
B22_FV=$(get_B22_stdout $outfile)

compare_floats $E_FV $refE_FV $E_tol
check_E_FV=$?
compare_floats $B20_FV $refB20_FV $B_tol
check_B20_FV=$?
compare_floats $B22_FV $refB22_FV $B_tol
check_B22_FV=$?

#-------------------------------------------------------------------------------
echo " --> Step 2: EFA calculation"

cat <<EOF >MOCCa.Cr55.data
&nucleus
neutrons=31, protons=24
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
blocklowest='n-'
/
&evolution
maxiter=$iterations 
/
&scfiteration
/
&wfs
nwn = 30, nwp = 30
/
&IO
InputFilename='MOCCa.Cr55.FV.wf'
Outputfilename='MOCCa.Cr55.EFA.wf'
/
&MomentParam
/
&Cranking
/
EOF

# Running the code
./$exebase <MOCCa.Cr55.data >../logs/blocking_EFA.out
mocca_check_EFA=$?
outfile="../logs/blocking_EFA.out"

# Check results for Step 2
E_EFA=$(get_total_energy_stdout $outfile)
B20_EFA=$(get_B20_stdout $outfile)
B22_EFA=$(get_B22_stdout $outfile)

compare_floats $E_EFA $refE_EFA $E_tol
check_E_EFA=$?
compare_floats $B20_EFA $refB20_EFA $B_tol
check_B20_EFA=$?
compare_floats $B22_EFA $refB22_EFA $B_tol
check_B22_EFA=$?

#-------------------------------------------------------------------------------
echo " --> Step 3: EFA blocking without T calculation"
exet="MOCCa.${exesuffix_T}.exe"

cat <<EOF >MOCCa.Cr55.data
&nucleus
neutrons=31, protons=24
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
blocklowest='n-'
/
&evolution
maxiter=$iterations 
/
&scfiteration
/
&wfs
nwn = 60, nwp = 60
/
&IO
InputFilename='MOCCa.Cr55.EFA.wf'
Outputfilename='MOCCa.Cr55.EFA.T.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF

# Running the code
./$exet <MOCCa.Cr55.data >../logs/blocking_EFA_T.out
mocca_check_EFA_T=$?
outfile="../logs/blocking_EFA_T.out"

# Check results for Step 3
E_EFA_T=$(get_total_energy_stdout $outfile)
B20_EFA_T=$(get_B20_stdout $outfile)
B22_EFA_T=$(get_B22_stdout $outfile)

compare_floats $E_EFA_T $refE_EFA_T $E_tol
check_E_EFA_T=$?
compare_floats $B20_EFA_T $refB20_EFA_T $B_tol
check_B20_EFA_T=$?
compare_floats $B22_EFA_T $refB22_EFA_T $B_tol
check_B22_EFA_T=$?

#-------------------------------------------------------------------------------
echo " --> Step 4: Full blocking"

cat <<EOF >MOCCa.Cr55.data
&nucleus
neutrons=31, protons=24
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
blocklowest='n-'
/
&evolution
maxiter=$iterations 
/
&scfiteration
/
&wfs
nwn = 60, nwp = 60
/
&IO
allowtransform=.true.
InputFilename='MOCCa.Cr55.EFA.wf'
Outputfilename='MOCCa.Cr55.block.wf'
/
&MomentParam
/
&Cranking
/
EOF

# Running the code
./$exet <MOCCa.Cr55.data >../logs/blocking_full.out
mocca_check_block=$?
outfile="../logs/blocking_full.out"

# Check results for Step 4
E_block=$(get_total_energy_stdout $outfile)
B20_block=$(get_B20_stdout $outfile)
B22_block=$(get_B22_stdout $outfile)

compare_floats $E_block $refE_block $E_tol
check_E_block=$?
compare_floats $B20_block $refB20_block $B_tol
check_B20_block=$?
compare_floats $B22_block $refB22_block $B_tol
check_B22_block=$?

#-------------------------------------------------------------------------------
# Clean up the work directory
cd ..
rm -r work/

#-------------------------------------------------------------------------------
# Print summary of all results
echo ""
echo "==============================================================================="
echo "                          BLOCKING TEST SUMMARY"
echo "==============================================================================="
echo ""

all_pass=true

# Step 1 summary
echo "Step 1: False Vacuum Calculation"
echo "  Expected Energy:  $refE_FV MeV"
echo "  Actual Energy:    $E_FV MeV"
echo "  Expected B20:     $refB20_FV"
echo "  Actual B20:       $B20_FV"
echo "  Expected B22:     $refB22_FV"
echo "  Actual B22:       $B22_FV"
if [ $check_E_FV -eq 0 ] && [ $check_B20_FV -eq 0 ] && [ $check_B22_FV -eq 0 ] && [ $mocca_check_FV -eq 0 ]; then
  echo "  Status:           PASS"
else
  echo "  Status:           FAIL"
  all_pass=false
fi
echo ""

# Step 2 summary
echo "Step 2: EFA Calculation"
echo "  Expected Energy:  $refE_EFA MeV"
echo "  Actual Energy:    $E_EFA MeV"
echo "  Expected B20:     $refB20_EFA"
echo "  Actual B20:       $B20_EFA"
echo "  Expected B22:     $refB22_EFA"
echo "  Actual B22:       $B22_EFA"
if [ $check_E_EFA -eq 0 ] && [ $check_B20_EFA -eq 0 ] && [ $check_B22_EFA -eq 0 ] && [ $mocca_check_EFA -eq 0 ]; then
  echo "  Status:           PASS"
else
  echo "  Status:           FAIL"
  all_pass=false
fi
echo ""

# Step 3 summary
echo "Step 3: EFA Blocking without T"
echo "  Expected Energy:  $refE_EFA_T MeV"
echo "  Actual Energy:    $E_EFA_T MeV"
echo "  Expected B20:     $refB20_EFA_T"
echo "  Actual B20:       $B20_EFA_T"
echo "  Expected B22:     $refB22_EFA_T"
echo "  Actual B22:       $B22_EFA_T"
if [ $check_E_EFA_T -eq 0 ] && [ $check_B20_EFA_T -eq 0 ] && [ $check_B22_EFA_T -eq 0 ] && [ $mocca_check_EFA_T -eq 0 ]; then
  echo "  Status:           PASS"
else
  echo "  Status:           FAIL"
  all_pass=false
fi
echo ""

# Step 4 summary
echo "Step 4: Full Blocking"
echo "  Expected Energy:  $refE_block MeV"
echo "  Actual Energy:    $E_block MeV"
echo "  Expected B20:     $refB20_block"
echo "  Actual B20:       $B20_block"
echo "  Expected B22:     $refB22_block"
echo "  Actual B22:       $B22_block"
if [ $check_E_block -eq 0 ] && [ $check_B20_block -eq 0 ] && [ $check_B22_block -eq 0 ] && [ $mocca_check_block -eq 0 ]; then
  echo "  Status:           PASS"
else
  echo "  Status:           FAIL"
  all_pass=false
fi
echo ""

# Overall summary
echo "==============================================================================="
if [ "$all_pass" = true ]; then
  echo "  OVERALL RESULT: ALL TESTS PASSED"
else
  echo "  OVERALL RESULT: SOME TESTS FAILED"
fi
echo "==============================================================================="
echo ""

# Return exit code 1 if any of the checks failed
if [ "$all_pass" = true ]; then
  exit 0
else
  exit 1
fi
