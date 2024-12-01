#-------------------------------------------------------------------------------
# Testing template
#
# Owner                : wouter.ryssens@ulb.be
# Complexity           : low
# Reference commit hash: commit 24470533c34a47af11223514bbecd500b212151a
#-------------------------------------------------------------------------------
# These are the hardcoded answers 
refE=-128.513125 # Total energy of O16 in MeV

# Basic starting point of all testing scripts
source ../functions.sh

# Set up and navigate to a work directory for
setup_test_env "default" "SLy4"

# Create runtime data
cat << EOF > tant.data
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='SLy4'
/
&pairing
/
&evolution
maxiter=1000
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='wtf'
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation
./$exe < tant.data > tant.out
# .... and immediately check it errors reported back by Tantalus
tantalus_error_codes $?

# Get the total energy from the STDOUT file
E=$(get_total_energy_stdout tant.out)
# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $E $refE 0.001

