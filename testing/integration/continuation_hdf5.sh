#-------------------------------------------------------------------------------
# Perform a constrained calculation writing an HDF5 file and then a second one
# to check the quality of the codes warmstarting capabilities
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - total energy                        <equal between both files> 1 keV
#  - convergence criteria satisfied      <both times>
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash continuation_HDF5.sh [EXESUFFIX] [PARAM]
#
# where EXESUFFIX specifies the  suffix of the executable to be used
#            Example: "BXL" for "Tantalus.BXL.exe".
# and PARAM specifies the parameterisation to be used.
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: commit 269c4b416f9682d4b0bef9fe321d9e24418f7a4
#-------------------------------------------------------------------------------

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

set -v
# Set up
setup_test_env "continuation_HDF5" "$1" "$2"

# Create runtime data for both code runs
cat << EOF > tant.data
&nucleus
neutrons=24, protons=24
energy_prec=1e-10
moment_prec=0.001
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$2"
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
osc_freq = 0.2, 0.2, 0.16
/
&IO
InputFilename='init'
OutputFilename='mf.hdf5'
/
&MomentParam
moreconstraints=.true.
/
&momentconstraint
l=2
m=0
constraint=50
moreconstraints=.true.
/
&momentconstraint
l=2
m=2
constraint=20
/
&Cranking
/
EOF

cat << EOF > tant.continuation.data
&nucleus
neutrons=24, protons=24
energy_prec=1e-07
moment_prec=0.001
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param="$2"
/
&pairing
type='HFB'
/
&evolution
maxiter=10
/
&scfiteration
/
&wfs
nwn = 50, nwp = 50
/
&IO
InputFilename='mf.hdf5'
OutputFilename='trash'
/
&MomentParam
moreconstraints=.true.
/
&momentconstraint
l=2
m=0
constraint=50
moreconstraints=.true.
multfromfile=.true.
/
&momentconstraint
l=2
m=2
constraint=20
multfromfile=.true.
/
&Cranking
/
EOF

# Run the first calculation from scratch
./$exe < tant.data > $outfile.a
tantalus_check=$?
# ... and then a second continuation calculation
./$exe < tant.continuation.data > $outfile.b
tantalus_check_continuation=$?

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Get the total energy from both files
Ea=$(get_total_energy_stdout $outfile.a)
Eb=$(get_total_energy_stdout $outfile.b)
# ... and compare with a tolerance of 1 keV to the expected answer
compare_floats $Ea $Eb 0.001
check_energy=$?
# ... and convergence!
check_convergence $outfile.b
c_b=$?
# ... and convergence!
check_convergence $outfile.a
c_a=$?

# remove working directory and traces of these calculations
teardown_test_env
#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
fail=$(($tantalus_check || $tantalus_check_continuation || $check_energy || $c_a || $c_b ))

if (($fail == 0)) ; then
    echo -e "test continuation_HDF5 :\033[1;32m success \033[0m"
else
    echo -e "test continuation_HDF5 :\033[1;31m failed ! 1st : $tantalus_check, 2nd : $tantalus_check_continuation, E : $check_E, convergence : $c_a $c_b \033[0m"
fi

exit $fail
