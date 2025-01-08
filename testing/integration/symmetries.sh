#---------------------------------------------------------------------------------
# Testing the equivalence of different symmetry modes with multiple calculations
# of Mg24 with BSkG3 constrained to a triaxial shape.
#
# More precisely, this testing script runs the exact same calculation with four
# different Tantalus configuration files:
#    - BXL   : Brussels functional(s) in maximally symmetric mode
#    - BXL-P : Brussels functional(s) with broken parity
#    - BXL-T : Brussels functional(s) with broken time-reversal
#    - BXL-TP: Brussels functional(s) with broken time-reversal AND broken parity
#
# The script does not check the results w.r.t. a hard-coded answer but rather
# verifies that the outcomes of all calculations are identical.
#
#  Quantity                              Answer                     Tolerance
#  --------                              ------                     ---------
#  Total energy                       < the first result obtained >   1 keV
#
#
#  TODO:
#  - add other observables
#
# Owner                : wouter.ryssens@ulb.be
# Complexity           : medium
# Reference commit hash: ec5d48ee9e546d795cc74d6160bf7812cd1d3b93
#-------------------------------------------------------------------------------
set -e
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Setting up the reference calculation: BXL
setup_test_env "symmetries" "BXL.gfortran.NUCLEI" "BSkG3"

cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='BSkG3'
/
&pairing
Type='HFB'
/
&evolution
maxiter=400
/
&scfiteration
/
&wfs
nwn = 24, nwp = 24
osc_freq = 0.2, 0.19, 0.18
/
&IO
InputFilename='init'
OutputFilename='trash'
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=20
moreconstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=10
/
&Cranking
/
EOF

# Run the calculation
./$exe < tant.data > $outfile
# .... and immediately check if Tantalus reported back some error codes
tantalus_check=$?
# Saving reference values
refE=$(get_total_energy_stdout $outfile)
# ... and tear down this testing environment.
teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Performing a parity-broken calculation
setup_test_env "symmetries" "BXL-P.gfortran.NUCLEI" "BSkG3"

cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
/
&mesh
nx=12, ny=12, nz=24, dx=1.0
/
&func
name_param='BSkG3'
/
&pairing
Type='HFB'
/
&evolution
maxiter=400
/
&scfiteration
/
&wfs
nwn = 24, nwp = 24
osc_freq = 0.2, 0.19, 0.18
/
&IO
InputFilename='init'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=20
moreconstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=10
/
&Cranking
/
EOF

# Run the calculation
./$exe < tant.data > $outfile
# ... immediately check if Tantalus reported back some error codes
tantalus_check_P=$?
# Saving reference values
E_parity=$(get_total_energy_stdout $outfile)
# ... remove all trace of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Performing a time-reversal-broken calculation
setup_test_env "symmetries" "BXL-T.gfortran.NUCLEI" "BSkG3"

cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
/
&mesh
nx=12, ny=12, nz=12, dx=1.0
/
&func
name_param='BSkG3'
/
&pairing
Type='HFB'
/
&evolution
maxiter=400
/
&scfiteration
/
&wfs
nwn = 48, nwp = 48
osc_freq = 0.2, 0.19, 0.18
/
&IO
InputFilename='init'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=20
moreconstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=10
/
&Cranking
/
EOF
# Run the calculation
./$exe < tant.data > $outfile
# ... immediately check if Tantalus reported back some error codes
tantalus_check_T=$?
# Saving reference values
E_timereversal=$(get_total_energy_stdout $outfile)
# ... remove all trace of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Performing a time-reversal- and parity-broken calculation
setup_test_env "symmetries" "BXL-P.gfortran.NUCLEI" "BSkG3"

# We need to set-up an initial run
cat << EOF > tant.init.data
&nucleus
neutrons=12, protons=12
/
&mesh
nx=12, ny=12, nz=24, dx=1.0
/
&func
name_param='BSkG3'
/
&pairing
Type='HFB'
/
&evolution
maxiter=0
/
&scfiteration
/
&wfs
nwn = 24, nwp = 24
osc_freq = 0.2, 0.19, 0.18
/
&IO
InputFilename='init'
OutputFilename='setting_up.wf'
allowtransform=.true.
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=20
moreconstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=10
/
&Cranking
/
EOF
# Run the initial run calculation (we don't keep the output)
./$exe < tant.init.data > /dev/null
# ... but we do keep the wf file!
mv setting_up.wf ../
teardown_test_env

# Change the environment
setup_test_env "symmetries" "BXL-TP.gfortran.NUCLEI" "BSkG3"
mv ../setting_up.wf .

# ... and now put the actual data!
cat << EOF > tant.data
&nucleus
neutrons=12, protons=12
/
&mesh
nx=12, ny=12, nz=24, dx=1.0
/
&func
name_param='BSkG3'
/
&pairing
Type='HFB'
/
&evolution
maxiter=400
/
&scfiteration
/
&wfs
nwn = 48, nwp = 48
osc_freq = 0.2, 0.19, 0.18
/
&IO
InputFilename='setting_up.wf'
OutputFilename='trash'
allowtransform=.true.
/
&MomentParam
moreconstraints=.true.
/
&MomentConstraint
l=2
m=0
constraint=20
moreconstraints=.true.
/
&MomentConstraint
l=2
m=2
constraint=10
/
&Cranking
/
EOF
# Run the calculation
./$exe < tant.data > $outfile
# ... immediately check if Tantalus reported back some error codes
tantalus_check_TP=$?
# Saving reference values
E_timereversal_parity=$(get_total_energy_stdout $outfile)
# ... remove all trace of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking
# a) Compare total energies with a tolerance of 1 keV
# a.1) parity-broken calculation
compare_floats $E_parity              $refE 0.001
check_energy_P=$?
# a.2) time-reversal broken calculation
compare_floats $E_timereversal        $refE 0.001
check_energy_T=$?
# a.3) time-reversal broken calculation
compare_floats $E_timereversal_parity $refE 0.001
check_energy_TP=$?

#- - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
t_check=$tantalus_check || $tantalus_check_P || $tantalus_check_T || $tantalus_check_TP
e_check=$check_energy_P || $check_energy_T || $check_energy_TP
return $t_check || $e_check
