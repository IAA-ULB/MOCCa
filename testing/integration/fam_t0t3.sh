#-------------------------------------------------------------------------------
# Perform a spherical HF + FAM calculation of O16 with t0t3 in a minimal box and 
# compare the Q_20 strength @ 25 MeV.
# This test may not immediately work for you, since it runs an executable called 
# Tantalus.LO.master.exe which is a LO compilation of the master branch, required 
# to be able to read in .wf file for the subsequent FAM calculation. 
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# note : the FAM code is sensitive to tiny changes in the reference state. Hence,
#        HF is converged up to high precision (E_prec = 1e-16). In addition, the
#        heavy-ball parameters (dt, mu) are fixed to ensure reproducability. Not
#        doing so would lead to unstable strengths differing from run to run.  
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total HF energy                -177.062001 MeV                  1 keV
#  - strength S_20 @ 25.0 MeV          1.665536 fm^4 MeV^-1          1e-6
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_t0t3.sh 
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : P. Demol [pepijn.demol@ulb.be]
# Reference commit hash: cd78f02f 
#-------------------------------------------------------------------------------
# These are the hardcoded answers
refE=-177.062001 # Total energy of O16 in MeV
refS20=1.665536  # Q_20 strength of O16 at 25 MeV in fm^4 MeV^-1  
# /!\ : this reference value is invalid as it differs on a run by run basis. 


# set -e # exit immediately if command gives non-zero exit status => disabled
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env_fam "fam_t0t3" "LO.master" "LO-T" "t0t3"

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
maxiter=1000
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
# (2) Run the FAM calculation

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
maxiter=10000
XY_prec=1e-10
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

# ... and compare with a tolerance of 1e-6 to the expected answer
compare_floats $S $refS20 0.000001
check_strength=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# remove working directory and traces of these calculations
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
fail=$(($tantalus_check || $fam_check || $check_energy || $check_strength ))

if (($fail == 0)) ; then
	echo -e "test FAM t0t3 :\033[1;32m success \033[0m"
else
	echo -e "test FAM t0t3 :\033[1;31m failed ! tant : $tantalus_check, fam : $fam_check, E_hf : $check_energy, S20 : $check_strength \033[0m"
fi

exit $fail
