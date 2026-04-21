#!/usr/bin/env sh
#--------------------------------------------------------------------------------
# Verify that the isoscalar quadrupole response obtained in (Q)FAM calculations
# satisfies
#    S(-\omega,+\eta) = S^*(+\omega,+\eta)
#    S(+\omega,-\eta) = S^*(+\omega,+\eta)
#    S(-\omega,-\eta) = S  (+\omega,+\eta)
#
# for O18 in a small modelspace.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
# bash fam_strenth_symmetry.sh [EXESUFFIX] [PARAM] [PAIRING] [OMEGA] [ETA] [-v/--verbose] [--free]
#
# where
#  - EXESUFFIX specifies the executable to be used
#  - PARAM the parameterisation.
#  - PAIRING = 'HF' or 'HFB'
#  - OMEGA   = the frequency [MeV]
#  - ETA     = the smearing  [MeV]
#
# Calling the script with flags
#    -v or --verbose will print the values which are compared and details on symmetry
#    --free will execute the test for the free response, i.e. without induced mean-fields
# Dependencies:
#   None
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner   : W. Ryssens [wouter.ryssens@ulb.be]
#--------------------------------------------------------------------------------
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh
LC_NUMERIC="en_US.UTF-8" # Important to deal with floats and particularly with scientific notation

# Initialize verbose mode as false by default
verbose=false
free=false

# Temporary array to hold arguments
args=()

# Parse all arguments for -v or --verbose
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--verbose)
      verbose=true
      shift
      ;;
    --free)
      free=true
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done


# Set up FAM testing environment
setup_test_env_fam "fam_strength_symmetry.${args[2]}.om=${args[3]}.eta=${args[4]}" "${args[0]}" "${args[0]}" "${args[1]}"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1a) Run a standard mean-field calculation
# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=10, protons=8
energy_prec=1e-16
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param="${args[1]}"
/
&pairing
type="${args[2]}"
/
&evolution
maxiter=100
dt=0.0209, momentum=0.5746
Estimateparams=.false.
/
&scfiteration
/
&wfs
nwn = 28, nwp = 28
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='mf_hfb.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation
./$exe < mf.data > $mfoutfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1b) Run a calculation that freezes the potentials just to get a
#     robust set of virtual states
cat << EOF > mf.data
&nucleus
neutrons=10, protons=8
energy_prec=1e-20
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param="${args[1]}"
/
&pairing
type="${args[2]}"
/
&evolution
maxiter=1000
dt=0.0209, momentum=0.5746
Estimateparams=.false.
freezeiter=1000
/
&scfiteration
/
&wfs
nwn = 28, nwp = 28
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='mf_hfb.wf'
OutputFilename='mf_hfb.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF
./$exe < mf.data > $mfoutfile.bis

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2c) Run the LO QFAM calculations

for smear in + -
do
for omega in + -
do

# Create runtime data
cat << EOF > qfam.data
&nucleus
neutrons=10, protons=8
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type="${args[2]}"
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 28, nwp = 28
/
&IO
InputFilename='mf_hfb.wf'
OutputFilename='trash'
famfile='S_20.QFAM.$omega.$smear.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=$omega${args[3]}
smear=$smear${args[4]}
l=2
m=0
maxiter=$(if [ "$free" = true ]; then echo "0"; else echo "30"; fi)
fam_precision=1e-8
/
EOF

#Redefine famoutfile to keep all logs
famoutfile=../logs/fam_strength_symmetry.${args[2]}.om=${args[3]}.eta=${args[4]}.${args[0]}.${args[1]}.fam.omega=$omega.$smear.out

# Run the calculation
./$exefam < qfam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
#qfam_check_eff1=$?

done
done

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Parse and compare the FAM output files
# Define tolerance for comparison
tolerance=0.0000001

# Extract strength components from all four files
components_pp=($(extract_strength_components "S_20.QFAM.+.+.fam"))
components_pm=($(extract_strength_components "S_20.QFAM.+.-.fam"))
components_mp=($(extract_strength_components "S_20.QFAM.-.+.fam"))
components_mm=($(extract_strength_components "S_20.QFAM.-.-.fam"))

# Assign components to meaningful variable names
# components_pp: S(ω, η) = S_re_pp + i*S_im_pp
# components_pm: S(ω, -η) = S_re_pm + i*S_im_pm
# components_mp: S(-ω, η) = S_re_mp + i*S_im_mp
# components_mm: S(-ω, -η) = S_re_mm + i*S_im_mm
S_re_pp=${components_pp[0]}
S_im_pp=${components_pp[1]}
S_re_pm=${components_pm[0]}
S_im_pm=${components_pm[1]}
S_re_mp=${components_mp[0]}
S_im_mp=${components_mp[1]}
S_re_mm=${components_mm[0]}
S_im_mm=${components_mm[1]}

# Print values if verbose mode is enabled
if [ "$verbose" = true ]; then
    echo "Strength function components:"
    printf "S( ω,  η) = %12.8f + i*%12.8f\n" $S_re_pp $S_im_pp
    printf "S( ω, -η) = %12.8f + i*%12.8f\n" $S_re_pm $S_im_pm
    printf "S(-ω,  η) = %12.8f + i*%12.8f\n" $S_re_mp $S_im_mp
    printf "S(-ω, -η) = %12.8f + i*%12.8f\n" $S_re_mm $S_im_mm
    echo ""
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3a) Test symmetry 1: S(ω) = S(-ω)* (complex conjugate)
# This means:
# S_re(ω, η) = S_re(-ω, η)
# S_im(ω, η) = -S_im(-ω, η)

if [ "$verbose" = true ]; then
    echo "Testing symmetry 1: S(ω, η) = S(-ω, η)*"
    printf "  Checking: S_re(ω, η) == +S_re(-ω, η) [%12.8f == %12.8f]\n" $S_re_pp $S_re_mp
    printf "  Checking: S_im(ω, η) == -S_im(-ω, η) [%12.8f == %12.8f]\n" $S_im_pp $(echo "-1 * $S_im_mp" | bc)
    echo ""
fi

# Test real parts equality
compare_floats "$S_re_pp" "$S_re_mp" "$tolerance"
sym1_re=$?

# Test imaginary parts equality (with sign change)
S_im_mp_neg=$(echo "-1 * $S_im_mp" | bc)
compare_floats "$S_im_pp" "$S_im_mp_neg" "$tolerance"
sym1_im=$?

# Overall symmetry 1 test
if [ $sym1_re -eq 0 ] && [ $sym1_im -eq 0 ]; then
    sym1_pass=true
    echo "✓ Symmetry 1 PASSED: S(ω, η) = S(-ω, η)*"
else
    sym1_pass=false
    echo "✗ Symmetry 1 FAILED: S(ω, η) ≠ S(-ω, η)*"
    if [ $sym1_re -ne 0 ]; then
        echo "  Real parts don't match: $S_re_pp ≠ $S_re_mp"
    fi
    if [ $sym1_im -ne 0 ]; then
        echo "  Imaginary parts don't match: $S_im_pp ≠ $(echo "-1 * $S_im_mp" | bc)"
    fi
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3b) Test symmetry 2: S(ω, η) = S(ω, -η)* (complex conjugate)
# This means:
# S_re(ω, η) = S_re(ω, -η)
# S_im(ω, η) = -S_im(ω, -η)

if [ "$verbose" = true ]; then
    echo "Testing symmetry 2: S(ω, η) = S(ω, -η)*"
    printf "  Checking: S_re(ω, η) == +S_re(ω, -η) [%12.8f == %12.8f]\n" $S_re_pp $S_re_pm
    printf "  Checking: S_im(ω, η) == -S_im(ω, -η) [%12.8f == %12.8f]\n" $S_im_pp $(echo "-1 * $S_im_pm" | bc)
    echo ""
fi

# Test real parts equality
compare_floats "$S_re_pp" "$S_re_pm" "$tolerance"
sym2_re=$?

# Test imaginary parts equality (with sign change)
S_im_pm_neg=$(echo "-1 * $S_im_pm" | bc)
compare_floats "$S_im_pp" "$S_im_pm_neg" "$tolerance"
sym2_im=$?

# Overall symmetry 2 test
if [ $sym2_re -eq 0 ] && [ $sym2_im -eq 0 ]; then
    sym2_pass=true
    echo "✓ Symmetry 2 PASSED: S(ω, η) = S(+ω,-η)*"
else
    sym2_pass=false
    echo "✗ Symmetry 2 FAILED: S(ω, η) ≠ S(+ω,-η)*"
    if [ $sym2_re -ne 0 ]; then
        echo "  Real parts don't match: $S_re_pp ≠ $S_re_pm"
    fi
    if [ $sym2_im -ne 0 ]; then
        echo "  Imaginary parts don't match: $S_im_pp ≠ $(echo "-1 * $S_im_pm" | bc)"
    fi
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3c) Test symmetry 3: S(ω, η) = S(-ω, -η) (no complex conjugate)
# This means:
# S_re(ω, η) = S_re(-ω, -η)
# S_im(ω, η) = S_im(-ω, -η)

if [ "$verbose" = true ]; then
    echo "Testing symmetry 3: S(ω, η) = S(-ω,-η)"
    printf "  Checking: S_re(ω, η) == S_re(-ω, -η) [%12.8f == %12.8f]\n" $S_re_pp $S_re_mm
    printf "  Checking: S_im(ω, η) == S_im(-ω, -η) [%12.8f == %12.8f]\n" $S_im_pp $S_im_mm
    echo ""
fi

# Test real parts equality
compare_floats "$S_re_pp" "$S_re_mm" "$tolerance"
sym3_re=$?

# Test imaginary parts equality
compare_floats "$S_im_pp" "$S_im_mm" "$tolerance"
sym3_im=$?

# Overall symmetry 3 test
if [ $sym3_re -eq 0 ] && [ $sym3_im -eq 0 ]; then
    sym3_pass=true
    echo "✓ Symmetry 3 PASSED: S(ω, η) = S(-ω,-η)"
else
    sym3_pass=false
    echo "✗ Symmetry 3 FAILED: S(ω, η) ≠ S(-ω,-η)"
    if [ $sym3_re -ne 0 ]; then
        echo "  Real parts don't match: $S_re_pp ≠ $S_re_mm"
    fi
    if [ $sym3_im -ne 0 ]; then
        echo "  Imaginary parts don't match: $S_im_pp ≠ $S_im_mm"
    fi
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Final summary
if [ "$sym1_pass" = true ] && [ "$sym2_pass" = true ] && [ "$sym3_pass" = true ]; then
    echo ""
    echo "=========================================="
    echo "ALL SYMMETRY TESTS PASSED ✓"
    echo "=========================================="
    exit 0
else
    echo ""
    echo "=========================================="
    echo "SOME SYMMETRY TESTS FAILED ✗"
    echo "=========================================="
    exit 1
fi
