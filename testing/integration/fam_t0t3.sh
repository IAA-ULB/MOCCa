#-------------------------------------------------------------------------------
# Perform a spherical HF + FAM calculation of O16 with t0t3 in a minimal box
#  - compare the Q_20 strength @ 25 MeV to a known result
#  - and redo it for an explicitly time-reversal broken FAM calculation
#  - and redo it for an explicitly parity broken FAM calculation
#  - and redo it for an explicitly particle-number broken QFAM calculation
#  - and redo it for an explicitly particle-number + time-reversal broken QFAM calculation
# In addition, you can specify a number of OpenMP threads to be used; if you
# do not provide an OpenMP-enabled executable, then this will simply be ignored.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total HF energy                -177.062001     MeV               1     keV
#  - strength S_20 @ 25.0 MeV          1.661155     fm^4 MeV^-1       1e-4  fm^4 MeV^-1
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# note : the FAM code is sensitive to tiny changes in the reference state. Hence,
#        HF is converged up to high precision (E_prec = 1e-16). In addition, the
#        heavy-ball parameters (dt, mu) are fixed to ensure reproducability. Not
#        doing so would lead to unstable strengths differing from run to run.
#
#        There are two different values the FAM code can converge to based on the
#        presence on the prediagonalisation of h at the start of fam_run.f90
#          - if the diagonalisation is enabled  : S_20(25.0) = 1.6611546044656
#          - if the diagonalisation is disabled : S_20(25.0) = 1.6612681654265
#         Note that the diagonalisation is done by default for HF calculations
#         and is implicit in HFB calculations with the two-basis method.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_t0t3.sh [EXESUFFIX] [EXESUFFIX_T] [EXESUFFIX_P] [Nthreads] [-v/--verbose]
#
# where EXESUFFIX,  EXESUFFIX_T, EXESUFFIX_P specify executables to be used:
# maximally symmetric, a time-reversal breaking and a parity breaking one.
# calling the script with flag -v or --verbose will print the values which are
# compared.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : P. Demol [pepijn.demol@ulb.be]
# Reference commit hash: bf93d26f7e11175189ffa93dad4d2065adf702eb
#-------------------------------------------------------------------------------

# Initialize verbose mode as false by default
verbose=false

# Temporary array to hold arguments
args=()

# Parse all arguments for -v or --verbose
while [[ $# -gt 0 ]]; do
  case "$1" in
    -v|--verbose)
      verbose=true
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

# These are the hardcoded answers
refE=-177.062001 # Total energy of O16 in MeV
refS20=1.661298165    # Q_20 strength of O16 at 25 MeV in fm^4 MeV^-1

if $verbose; then
	echo "reference values and tolerances: "
    echo "               E = $refE (+/- 0.001) MeV"
    echo "     S(omega=25) = $refS20 (+/- 0.001) fm^4 MeV^-1 "
    echo
fi

Nthreads=${args[3]}
export OMP_NUM_THREADS=$Nthreads
echo "Running with " $Nthreads " OpenMP threads"
set -e # exit immediately if command gives non-zero exit status
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env_fam "fam_t0t3_Nthreads=$Nthreads" "${args[0]}" "${args[0]}" "t0t3"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1) Run the mean-field calculation

# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=8, protons=8
energy_prec=1e-16
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type='HF'
/
&evolution
maxiter=100
dt=0.0209, momentum=0.5746
Estimateparams=.false.
/
&scfiteration
/
&wfs
nwn = 14, nwp = 14
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='mf.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation
./$exe < mf.data > $mfoutfile
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2) Run a calculation that freezes the potentials just to get a
#     robust set of virtual states
cat << EOF > mf.data
&nucleus
neutrons=8, protons=8
energy_prec=1e-20
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type='HF'
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
nwn = 14, nwp = 14
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='mf.wf'
OutputFilename='mf.wf'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
EOF
./$exe < mf.data > $mfoutfile.bis


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Run the LO FAM calculation

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 14, nwp = 14
/
&IO
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='S_20.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega_min=25
omega_max=25.0
smear=1.0
l=2
m=0
maxiter=30
fam_precision=1e-8
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if MOCCa reported back some error codes
fam_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the total energy from the STDOUT file
E=$(get_total_energy_stdout $mfoutfile)

if $verbose; then
	echo " HF  LO   : E = $E"
fi

# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $E $refE 0.001
check_energy=$?

# b) Get the strength from the S_20.fam file
S=$(get_strength "S_20.fam" 25.0)

if $verbose; then
	echo " FAM LO   : S(omega=25) = $S"
fi

# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S $refS20 0.001
check_strength=$?

# We keep the wf file!
mv mf.wf ../
# ... but otherwise clean-up
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Run the LO-T FAM calculation
setup_test_env_fam "fam_t0t3_Nthreads=$Nthreads" "${args[0]}" "${args[1]}" "t0t3"
cp ../mf.wf .

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
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
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='S_20.T.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=25.0
smear=1.0
l=2
m=0
maxiter=30
fam_precision=1e-8
/
EOF
# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if MOCCa reported back some error codes
fam_T_check=$?


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# b) Get the strength from the S_20.T.fam file
S_T=$(get_strength "S_20.T.fam" 25.0)

if $verbose; then
	echo " FAM LO-T : S(omega=25) = $S_T"
fi

# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S_T $refS20 0.001
check_strength_T=$?
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# remove working directory and traces of these calculations
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (4) Run the LO-P FAM calculation
setup_test_env_fam "fam_t0t3_Nthreads=$Nthreads" "${args[0]}" "${args[2]}" "t0t3"
mv ../mf.wf .

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=8, ny=8, nz=16, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type='HF'
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 14, nwp = 14
/
&IO
InputFilename='mf.wf'
OutputFilename='trash'
allowtransform=.true.
famfile='S_20.P.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=25.0
smear=1.0
l=2
m=0
maxiter=30
fam_precision=1e-10
/
EOF
# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if MOCCa reported back some error codes
fam_P_check=$?
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# b) Get the strength from the S_20.fam file
S_P=$(get_strength "S_20.P.fam" 25.0)

if $verbose; then
	echo " FAM LO-P : S(omega=25) = $S_P"
fi

# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S_P $refS20 0.001
check_strength_P=$?
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# remove working directory and traces of these calculations
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (5) Run the mean-field calculation with pairing (eventhough zero)
setup_test_env_fam "qfam_t0t3" "${args[0]}" "${args[0]}" "t0t3"

# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=8, protons=8
energy_prec=1e-16
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type='HFB'
/
&evolution
maxiter=100
dt=0.0209, momentum=0.5746
Estimateparams=.false.
/
&scfiteration
/
&wfs
nwn = 14, nwp = 14
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
./$exe < mf.data > $mfoutfile.ter
# .... and immediately check if MOCCa reported back some error codes
mocca_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (6) Run a calculation that freezes the potentials just to get a
#     robust set of virtual states
cat << EOF > mf.data
&nucleus
neutrons=8, protons=8
energy_prec=1e-20
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type='HFB'
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
nwn = 14, nwp = 14
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

E=$(get_total_energy_stdout $mfoutfile.ter)

if $verbose; then
	echo " HFB LO   : E = $E"
fi

# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $E $refE 0.001
check_energy=$?


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (7) Run the LO QFAM calculation

# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type='HFB'
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 14, nwp = 14
/
&IO
InputFilename='mf_hfb.wf'
OutputFilename='trash'
famfile='S_20.QFAM.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=25.0
smear=1.0
l=2
m=0
maxiter=30
fam_precision=1e-8
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if MOCCa reported back some error codes
qfam_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking

# b) Get the strength from the S_20.fam file
S_Q=$(get_strength "S_20.QFAM.fam" 25.0)

if $verbose; then
	echo " QFAM LO   : S(omega=25) = $S_Q"
fi

# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S_Q $refS20 0.001
check_strength_QRPA=$?

mv mf_hfb.wf ../


# ... but otherwise clean-up
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (8) Run the LO-T QFAM calculation

setup_test_env_fam "qfam_t0t3" "${args[0]}" "${args[1]}" "t0t3"

cp ../mf_hfb.wf .


# Create runtime data
cat << EOF > fam.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=8, ny=8, nz=8, dx=0.8
/
&func
name_param='t0t3'
/
&pairing
type='HFB'
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
famfile='S_20.QFAM.T.fam'
allowtransform=.true.
/
&MomentParam
/
&Cranking
/
&fam
omega=25.0
smear=1.0
l=2
m=0
maxiter=30
fam_precision=1e-8
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if MOCCa reported back some error codes
qfam_T_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking

# b) Get the strength from the S_20.fam file
S_QT=$(get_strength "S_20.QFAM.T.fam" 25.0)

if $verbose; then
	echo " QFAM LO-T : S(omega=25) = $S_QT"
fi


# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S_QT $refS20 0.001
check_strength_QRPA_T=$?

# ... but otherwise clean-up
teardown_test_env


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
fail=$(($mocca_check || $fam_check || $fam_T_check || $fam_P_check || $qfam_check || $qfam_T_check || $check_energy || $check_strength || $check_strength_T || $check_strength_P || $check_strength_QRPA || $check_strength_QRPA_T))

if (($fail == 0)) ; then
	echo -e "test FAM t0t3 :\033[1;32m success \033[0m"
	echo -e "  tantalus : $tantalus_check, fam    : $fam_check, fam_T : $fam_T_check, fam_P : $fam_P_check"
	echo -e "  E_hf     : $check_energy, S20    : $check_strength, S20_T : $check_strength_T, S20_P : $check_strength_P"
else
	echo -e "test FAM t0t3 :\033[1;31m failed ! exit status : tant = $mocca_check, fam = $fam_check, fam_T = $fam_T_check, fam_P = $fam_P_check, fam_HFB = $qfam_check
	                 benchmarks :  E_hf = $check_energy, S20 = $check_strength, S20_T = $check_strength_T, S20_P = $check_strength_P, S20_QRPA = $check_strength_QRPA, , S20_QRPA_T = $check_strength_QRPA_T \033[0m"
fi

exit $fail
