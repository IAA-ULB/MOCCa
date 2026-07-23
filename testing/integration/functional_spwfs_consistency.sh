#------------------------------------------------------------------------------------
# Test that the difference between energy calculated through the functional
# and by way of the spwfs (E_fu - E_sp) can be made smaller than 1e-11.
# This test uses the same parameters as run_mf.sh.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - |E_fu - E_sp|                        < 1e-11                    1e-11
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Usage
# -----
#   bash functional_spwfs_consistency.sh -e EXESUFFIX -p PARAM
#
# where:
#   -e EXESUFFIX: specifies the suffix of the executable to be used
#                 Example: "NLO" for "MOCCa.NLO.exe".
#   -p PARAM:    specifies the parameterization name (without .param extension)
#                 Example: "SLy4_NC"
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: (to be filled)
#------------------------------------------------------------------------------------

usage() { echo "Usage: $0 -e EXESUFFIX -p PARAM" 1>&2; exit 1; }

while getopts ":e:p:" opt; do
    case "${opt}" in
        e)
            exec_suffix=${OPTARG}
            ;;
        p)
            param=${OPTARG}
            ;;
        *)
            usage
            ;;
    esac
done

# Check that required arguments are provided
if [ -z "${exec_suffix}" ] || [ -z "${param}" ]; then
    usage
fi

set -e
#- - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

echo "============================================================"
echo "  Functional vs Single-Particle Energy Consistency Test"
echo "============================================================"
echo ""
echo "Configuration:"
echo "  Executable suffix : ${exec_suffix}"
echo "  Parameterization  : ${param}"
echo ""

# Set up
setup_test_env "functional_spwfs_consistency" "${exec_suffix}" "${param}"

echo "Running MOCCa calculation..."

# Create runtime data matching run_mf.sh
cat << EOF > mocca.data
&nucleus
neutrons=6
protons =6
energy_prec=1e-16
/
&mesh
nx=12, ny=12, nz=12
dx=0.8
/
&func
name_param="${param}"
/
&pairing
Type       ='HF'
/
&indices
/
&evolution
maxiter=1500
printiter=100
/
&scfiteration
/
&wfs
nwn =15, nwp =15
/
&IO
InputFilename='INIT'
Outputfilename='trash'
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation
./$exe < mocca.data > $outfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
echo ""
echo "Checking results..."
echo "-------------------"
echo ""

# a) Get the E_fu - E_sp difference from the STDOUT file
E_fu_diff=$(get_E_fu_diff_stdout $outfile)
echo "E_fu - E_sp difference : $E_fu_diff"

# Take absolute value and convert to decimal with enough precision in one awk command
# Use LC_NUMERIC to ensure decimal point is correctly interpreted
abs_E_fu_diff_decimal=$(LC_NUMERIC="en_US.UTF-8" echo "$E_fu_diff" | LC_NUMERIC="en_US.UTF-8" awk '{val = ($1 < 0 ? -$1 : $1); printf "%.20f", val}')
echo "Absolute difference    : $abs_E_fu_diff_decimal"

# b) Compare with a tolerance of 1e-11
compare_floats $abs_E_fu_diff_decimal 0 0.00000000001
check_E_fu_diff=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
echo ""
# Combine all checks
overall_check=$(($mocca_check || $check_E_fu_diff))

# Print final status
if [ $overall_check -eq 0 ]; then
    echo "All checks passed: SUCCESS"
else
    echo "Some checks failed: FAIL"
    echo ""
    echo "Detailed diagnostics:"
    echo "--------------------"
    
    # Specify which ones failed with detailed diagnostics
    if [ $mocca_check -ne 0 ]; then
        echo "  - MOCCa check failed: MOCCa returned non-zero exit code ($mocca_check)"
    fi
    
    if [ $check_E_fu_diff -ne 0 ]; then
        echo "  - Functional vs SPWFS consistency check failed:"
        echo "    Expected: < 1e-11"
        echo "    Obtained: $abs_E_fu_diff_decimal"
        echo "    Tolerance: 1e-11"
    fi
    
    echo ""
    echo "For more details, check the output file: $outfile"
fi

# Exit with the combined status
exit $overall_check
