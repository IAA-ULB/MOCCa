#-------------------------------------------------------------------------------
# Test the magnetic dipole moment of the 7/2- ground-state in Sc41
#
# This script tests the magnetic moment calculation for Sc41 using the SLy5s1
# parameterization in a two-step process:
#  1. False vacuum calculation (time-reversal conserving)
#  2. Blocked calculation with time-reversal breaking
#
# and compares the z-component of the magnetic moment (mu_z) to a known reference value.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - z-component of magnetic moment     5.80                       0.01
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash magmoment_Sc41.sh exec exec_T
#
# where
# - exec specifies the suffix of the executable to be used when T is conserved
#            Example: "NLO" for "MOCCa.NLO.exe".
#
# - exec_T specifies the suffix of the executable to be used when T is broken
#            Example: "NLO-T" for "MOCCa.NLO-T.exe".
#
# Dependencies: SLy5s1 parameterization
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit     : HEAD
#-------------------------------------------------------------------------------
# KNOWN RESULTS TO COMPARE AGAINST
ref_muz=5.80 # z-component of magnetic moment

# Tolerance for comparison
mu_tol=0.01

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

echo "==============================================================================="
echo "  Starting magnetic moment test for Sc41"
echo "==============================================================================="
echo ""

# Set up a single work directory for both steps
exesuffix=$1
exesuffix_T=$2
exebase="MOCCa.$exesuffix.exe"

echo "  Setting up test environment..."
if [ ! -d "work/" ]; then
  mkdir work
fi
if [ ! -d "logs/" ]; then
  mkdir logs
fi
if [ ! -d "../logs/" ]; then
  mkdir ../logs
fi

echo "  Copying executables and parameterization..."
cp ../../exec/$exebase work/
cp ../../exec/MOCCa.${exesuffix_T}.exe work/ 2>/dev/null || true
cp ../../parameterizations/SLy5s1.param work/

cd work

# Clean up any previous wavefunctions
rm -f MOCCa.Sc41.*.wf

echo ""
echo " --> Step 1: False vacuum calculation"

cat <<EOF >MOCCa.Sc41.data
&nucleus
neutrons=20, protons=21
/
&mesh
nx=12, ny=12, nz=12, dx=0.8
/
&func
name_param="SLy5s1"
/
&pairing
type='HFB'
/
&evolution
maxiter=100
printiter=10
/
&scfiteration
/
&wfs
nwn = 30, nwp = 30
/
&IO
InputFilename='init'
Outputfilename='MOCCa.Sc41.FV.wf'
/
&MomentParam
/
&MomentConstraint
/
&Cranking
/
EOF

# Running the code
./$exebase <MOCCa.Sc41.data >../logs/magmoment_Sc41_FV.out
mocca_check_FV=$?
outfile="../logs/magmoment_Sc41_FV.out"

if [ $mocca_check_FV -eq 0 ]; then
  echo "  Step 1 completed successfully"
else
  echo "  Step 1 failed with exit code $mocca_check_FV"
fi

#-------------------------------------------------------------------------------
echo ""
echo " --> Step 2: Blocked calculation"
exet="MOCCa.${exesuffix_T}.exe"

cat <<EOF >MOCCa.Sc41.data
&nucleus
neutrons=20, protons=21
/
&mesh
nx=12, ny=12, nz=12, dx=0.8
/
&func
name_param="SLy5s1"
/
&pairing
type='HFB'
blocknumber=1
blocktype=1    
/
&indices
blockindices=114 ! hard-coded blocking index
/
&evolution
maxiter=100
printiter=10
/
&scfiteration
/
&wfs
nwn = 60, nwp = 60
/
&IO
InputFilename='MOCCa.Sc41.FV.wf'
Outputfilename='MOCCa.Sc41.blocked.wf'
allowtransform=.true.
/
&MomentParam
/
&MomentConstraint
/
&Cranking
/
EOF

# Running the code
./$exet <MOCCa.Sc41.data >../logs/magmoment_Sc41_blocked.out
mocca_check_block=$?
outfile="../logs/magmoment_Sc41_blocked.out"

if [ $mocca_check_block -eq 0 ]; then
  echo "  Step 2 completed successfully"
else
  echo "  Step 2 failed with exit code $mocca_check_block"
fi

echo ""
echo "  Extracting magnetic moment from output..."

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the z-component of the magnetic moment from the STDOUT file
mu_z=$(get_muz_stdout $outfile)
# ... and compare with the tolerance to the expected answer
compare_floats $mu_z $ref_muz $mu_tol
check_muz=$?

#-------------------------------------------------------------------------------
# Clean up the work directory
cd ..
rm -r work/

#-------------------------------------------------------------------------------
# Print summary of the results
echo ""
echo "==============================================================================="
echo "                          MAGNETIC MOMENT TEST SUMMARY"
echo "==============================================================================="
echo ""
echo "Test: Magnetic dipole moment of 7/2- ground-state in Sc41"
echo "Parameterization: SLy5s1"
echo ""
echo "Expected mu_z:  $ref_muz"
echo "Actual mu_z:    $mu_z"
echo "Tolerance:      $mu_tol"
echo ""

all_pass=true

if [ $mocca_check_FV -eq 0 ] && [ $mocca_check_block -eq 0 ] && [ $check_muz -eq 0 ]; then
  echo "Status: PASS"
else
  echo "Status: FAIL"
  all_pass=false
fi

echo ""
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
