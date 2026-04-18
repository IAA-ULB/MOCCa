#!/usr/bin/env sh
#--------------------------------------------------------------------------------
# Verify that the strength obtained in FAM calculations satisfies multiple
#  symmetries.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
# bash fam_strenth_symmetry.sh [EXESUFFIX] [PARAM] [-v/--verbose]
#
# where EXESUFFIX specifies the executable to be used and PARAM the parameterisation.
# Calling the script with flag -v or --verbose will print the values which are
# compared.
#
#
# Dependencies:
#   None
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner   : W. Ryssens [wouter.ryssens@ulb.be]
#--------------------------------------------------------------------------------
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

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


# Set up FAM testing environment
setup_test_env_fam "fam_strength_symmetry" "${args[0]}" "${args[0]}" "${args[1]}"

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

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2c) Run the LO QFAM calculations

for smear in +1 -1
do
for omega in +25.0 -25.0
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
famfile='S_20.QFAM.$omega.$smear.fam'
/
&MomentParam
/
&Cranking
/
&fam
omega=$omega
smear=$smear
l=2
m=0
maxiter=30
fam_precision=1e-8
/
EOF

#Redefine famoutfile to keep all logs
famoutfile=fam_strength_symmetry.${args[0]}.${args[1]}.fam.omega=$omega.$smear.out

echo "Running FAM with output file " $famoutfile

# Run the calculation
./$exefam < qfam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
#qfam_check_eff1=$?

done
done
