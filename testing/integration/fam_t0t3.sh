#-------------------------------------------------------------------------------
# Perform a spherical HF + FAM calculation of O16 with t0t3 in a minimal box and 
# compare the Q_20 strength @ 25 MeV
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total HF energy                -177.062001 MeV                  1 keV
#  - strength S_20 @ 25.0MeV           1.632650 fm^4 MeV^-1          0.001
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_simple.sh 
#
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : P. Demol [pepijn.demol@ulb.be]
# Reference commit hash: fcbfa3c08956ff621ee58b6256a3cdd0ec14f270 
#-------------------------------------------------------------------------------
# These are the hardcoded answers
refE=-177.062001 # Total energy of O16 in MeV
refS20=1.632650
omega=25.000

set -e
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env_fam "fam_simple" "LO.master" "LO-T" "t0t3"

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# (1) Run the MF calculation

# Create runtime data
cat << EOF > mf.data
&nucleus
neutrons=8, protons=8
energy_prec=1e-12
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
nwn = 12, nwp = 12
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

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
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
nwn = 24, nwp = 24
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
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_check=$?
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the total energy from the STDOUT file
E=$(get_total_energy_stdout $mfoutfile)
# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $E $refE 0.001
check_energy=$?
# b) Get the strength from the S_20.fam file
S=$(get_strength S_20.fam $omega)
# ... and compare with a tolerance of 0.001 to the expected answer
compare_floats $S $refS20 0.001
check_strength=$?

# remove working directory and traces of these calculations
teardown_test_env
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
exit $(($tantalus_check || $fam_check || $check_energy || $check_strength ))
