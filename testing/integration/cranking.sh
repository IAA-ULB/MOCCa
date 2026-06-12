#-------------------------------------------------------------------------------
# Perform three cranked calculations for Cr48 with SLy4 
#  - one with a fixed Jz=2 hbar and direct diagonalisation for HFB
#  - one with a fixed omega_z=0.3292333881 MeV/hbar and direct diagonalisation
#  - one with a fixed Jz=2 hbar and gradient solver for HFB.
# All should give the same total energy of -411.090 MeV.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total energy                        -411.090 MeV                1 keV
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash cranking.sg [EXESUFFIX]
#
# where EXESUFFIX specifies the  suffix of the executable to be used
#            Example: "NLO-T" for "MOCCa.NLO-T.exe".
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: commit dd67b1bda68595904c5701032d95802e1236f546
#-------------------------------------------------------------------------------
# These are the hardcoded answers
refE=-411.090

set -e
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env "cranking" "$1" "SLy4"

# Create runtime data
cat << EOF > tant.data
&nucleus
neutrons=24, protons=24
energy_prec=1e-8
angmom_prec=1e-5
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy4'
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
nwn = 50, nwp = 50
osc_freq = 0.19, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
/
&Cranking
cranktypeZ=1
crankZ=2
/
EOF

# Run the calculation
./$exe < tant.data > $outfile.J
# .... and immediately check if MOCCa reported back some error codes
mocca_check_J=$?

# Second run: constant omega_z and direct diagonalisation
cat << EOF > tant.data
&nucleus
neutrons=24, protons=24
energy_prec=1e-8
angmom_prec=1e-5
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy4'
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
nwn = 50, nwp = 50
osc_freq = 0.19, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
/
&Cranking
omegaZ=0.3292333881
/
EOF
./$exe < tant.data > $outfile.omega
mocca_check_omega=$?

# Create runtime data
cat << EOF > tant.data
&nucleus
neutrons=24, protons=24
energy_prec=1e-8
angmom_prec=1e-5
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy4'
/
&pairing
type='HFB'
pairingscheme=1
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 50, nwp = 50
osc_freq = 0.19, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
/
&Cranking
cranktypeZ=1
crankZ=2
/
EOF

# Run the calculation
./$exe < tant.data > $outfile.gradient
# .... and immediately check if MOCCa reported back some error codes
mocca_check_gradient=$?

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the total energy from the STDOUT file
E_J=$(get_total_energy_stdout $outfile.J)
E_omega=$(get_total_energy_stdout $outfile.omega)
E_gradient=$(get_total_energy_stdout $outfile.gradient)
# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $E_J $refE 0.001
check_energy_J=$?
compare_floats $E_omega $refE 0.001
check_energy_omega=$?
compare_floats $E_gradient $refE 0.001
check_energy_gradient=$?
# remove working directory and traces of these calculations
teardown_test_env
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
exit $(( $mocca_check_J || $check_energy_J || $mocca_check_omega || $check_energy_omega || $mocca_check_gradient || $check_energy_gradient ))
