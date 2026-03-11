#-------------------------------------------------------------------------------
# Verify the linearity of the (Q)FAM by repeating a calculation with effective 
# charges 1.0 and 2.0 . This factor 2.0 should lead to a factor 4.0 in the total 
# strength. This is tested for: 
#  - FAM : O16 t0t3 in a minimal box, Q_20 strength @ 25 MeV
#  - QFAM : O18 t0t3 in a mimimal box, Q_20 strength @ 25 MeV
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash fam_t0t3.sh
#
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : P. Demol [pepijn.demol@ulb.be]
#-------------------------------------------------------------------------------

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env_fam "fam_t0t3" "LO" "LO" "t0t3"

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1a) Run the mean-field calculation

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
# (1b) Run a calculation that freezes the potentials just to get a 
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
# (1c) Run the LO FAM calculation

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
famfile='S_20.effch1.fam'
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
eff_charge_n=1.0
eff_charge_p=1.0
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_check_eff1=$?


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (1d) Repeat LO FAM calculation with effective charges = 2.0

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
famfile='S_20.effch2.fam'
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
eff_charge_n=2.0
eff_charge_p=2.0
/
EOF

# Run the calculation
./$exefam < fam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
fam_check_eff2=$?


# Get the strength from the S_20.fam file
Seff1=$(get_strength "S_20.effch1.fam" 25.0)
Seff2=$(get_strength "S_20.effch2.fam" 25.0)

# devide the strength of the last results by four
Seff2_scaled=$(echo "0.25 * $Seff2" | bc)



# ... and compare with a tolerance of 1e-4 to the expected answer
compare_floats $Seff1 $Seff2_scaled 0.0001
check_strength=$?

teardown_test_env

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2a) Run the HFB calculation for O18 
setup_test_env_fam "qfam_t0t3" "LO" "LO" "t0t3"

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
# .... and immediately check if Tantalus reported back some error codes
tantalus_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2b) Run a calculation that freezes the potentials just to get a 
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


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2c) Run the LO QFAM calculation

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
famfile='S_20.QFAM.effch1.fam'
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
eff_charge_n=1.0
eff_charge_p=1.0
/
EOF

# Run the calculation
./$exefam < qfam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
qfam_check_eff1=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# (2d) Repeat LO QFAM calculation with effective charges 2.0

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
famfile='S_20.QFAM.effch2.fam'
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
eff_charge_n=2.0
eff_charge_p=2.0
/
EOF

# Run the calculation
./$exefam < qfam.data > $famoutfile
# .... and immediately check if Tantalus reported back some error codes
qfam_check_eff2=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Starting the checking

# Get the strength from the S_20.fam file
Seff1=$(get_strength "S_20.QFAM.effch1.fam" 25.0)
Seff2=$(get_strength "S_20.QFAM.effch2.fam" 25.0)

# devide the strength of the last results by four
Seff2_scaled=$(echo "0.25 * $Seff2" | bc)


# ... and compare with a tolerance of 1e-4 to the expected answer
compare_floats $Seff1 $Seff2_scaled 0.0001
check_strength_QRPA=$?

# ... but otherwise clean-up
teardown_test_env


#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
fail=$(($tantalus_check || $fam_check_eff1 || $fam_check_eff2 || $qfam_check_eff1 || $qfam_check_eff2 || $check_strength || $check_strength_QRPA ))

if (($fail == 0)) ; then
	echo -e "test FAM linearity :\033[1;32m success \033[0m"
else
	echo -e "test FAM linearity :\033[1;31m failed ! exit status : tant = $tantalus_check, fam_eff1 = $fam_check_eff1, fam_eff2 = $fam_check_eff1, qfam_eff1 = $qfam_check_eff1, qfam_eff2 = $qfam_check_eff1,
	                 benchmarks :  FAM S20(eff=1) == S20(eff=2) / 4 : $check_strength
	               	  QFAM S20(eff=1) == S20(eff=2) / 4 : $check_strength_QRPA \033[0m"
fi

exit $fail
