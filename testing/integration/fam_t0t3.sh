#-------------------------------------------------------------------------------
# Perform a spherical HF + FAM calculation of O16 with t0t3 in a minimal box
#  - compare the Q_20 strength @ 25 MeV to a known result
#  - and redo it for an explicitly time-reversal broken FAM calculation
#  - and redo it for an explicitly parity broken FAM calculation
# 
# In addition, you can specify a number of OpenMP threads to be used; if you
# do not provide an OpenMP-enabled executable, then this will simply be ignored.
# 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# note : the FAM code is sensitive to tiny changes in the reference state. Hence,
#        HF is converged up to high precision (E_prec = 1e-16). In addition, the
#        heavy-ball parameters (dt, mu) are fixed to ensure reproducability. Not
#        doing so would lead to unstable strengths differing from run to run.  
#
#        It turns out this is not even entirely sufficient; this is why the
#        test accuracy on the strength has been drastically reduced.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total HF energy                -177.062001 MeV                  1 keV
#  - strength S_20 @ 25.0 MeV          1.66     fm^4 MeV^-1          1e-2
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_t0t3.sh [EXESUFFIX] [EXESUFFIX_T] [EXESUFFIX_P] [Nthreads]
#
# where EXESUFFIX and EXESUFFIX_T specify executables to be used: a time-reversal
# conserving and a time-reversal breaking one.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : P. Demol [pepijn.demol@ulb.be]
# Reference commit hash: bf93d26f7e11175189ffa93dad4d2065adf702eb
#-------------------------------------------------------------------------------
# These are the hardcoded answers
refE=-177.062001 # Total energy of O16 in MeV
refS20=1.66      # Q_20 strength of O16 at 25 MeV in fm^4 MeV^-1

Nthreads=$4
export OMP_NUM_THREADS=$Nthreads
set -e # exit immediately if command gives non-zero exit status
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env_fam "fam_t0t3_Nthreads=$Nthreads" "$1" "$1" "t0t3"

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
# .... and immediately check if Tantalus reported back some error codes
tantalus_check=$?
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
# .... and immediately check if Tantalus reported back some error codes
fam_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the total energy from the STDOUT file
E=$(get_total_energy_stdout $mfoutfile)
# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $E $refE 0.001
check_energy=$?

# b) Get the strength from the S_20.fam file
S=$(get_strength "S_20.fam" 25.0)
# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S $refS20 0.01
check_strength=$?

# We keep the wf file!
mv mf.wf ../
# ... but otherwise clean-up
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (3) Run the LO-T FAM calculation
setup_test_env_fam "fam_t0t3_Nthreads=$Nthreads" "$1" "$2" "t0t3"
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
famfile='S_20.fam'
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
# .... and immediately check if Tantalus reported back some error codes
fam_T_check=$?


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# b) Get the strength from the S_20.fam file
S_T=$(get_strength "S_20.fam" 25.0)
# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S_T $refS20 0.01
check_strength_T=$?
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# remove working directory and traces of these calculations
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (4) Run the LO-P FAM calculation
setup_test_env_fam "fam_t0t3_Nthreads=$Nthreads" "$1" "$3" "t0t3"
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
# .... and immediately check if Tantalus reported back some error codes
fam_P_check=$?
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# b) Get the strength from the S_20.fam file
S_P=$(get_strength "S_20.fam" 25.0)
# ... and compare with a tolerance of 1e-2 to the expected answer
compare_floats $S_P $refS20 0.01
check_strength_P=$?
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# remove working directory and traces of these calculations
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
fail=$(($tantalus_check || $fam_check || $fam_T_check || $fam_P_check || $check_energy || $check_strength || $check_strength_T || $check_strength_P))

if (($fail == 0)) ; then
	echo -e "test FAM t0t3 :\033[1;32m success \033[0m"
	echo -e "  tantalus : $tantalus_check, fam    : $fam_check, fam_T : $fam_T_check"
	echo -e "  E_hf     : $check_energy, S20    : $check_strength, S20_T : $check_strength_T"
else
	echo -e "test FAM t0t3 :\033[1;31m failed ! tant : $tantalus_check, fam : $fam_check, E_hf : $check_energy, S20 : $check_strength, S20_T : $check_strength_T, S20_P : $check_strength_P \033[0m"
fi

exit $fail
