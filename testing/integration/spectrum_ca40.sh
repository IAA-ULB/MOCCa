#-------------------------------------------------------------------------------
# Frozen-potential spectrum calculation of Ca40 with SLy4, Hartree-Fock
# (PairingType 0), starting from a converged HF wavefunction file and
# solving the discretized spectrum up to Fermi + 20 MeV excitation.
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  Quantity                              Target                     Tolerance
#  --------                              ------                     ---------
#  - Chebyshev sizing runs               SPECTRUM SIZING banner     present
#  - certified straddle check            no error at end of run     pass
#  - trim to the cutoff                  Trimming output spectrum   present
#  - total energy (trimmed Ritz set)     = HF energy of run 1       1 keV
#
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash spectrum_ca40.sh [EXESUFFIX]
#
# where EXESUFFIX specifies the  suffix of the executable to be used
#            Example: "BXL" for "MOCCa.BXL.exe".
#
# Dependencies: none
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: n/a (item 007: Chebyshev count + fixed-block LOBPCG)
#-------------------------------------------------------------------------------

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# ENSURE small number of THREADS
export OMP_NUM_THREADS=4

# Set up
setup_test_env "spectrum_ca40" "$1" "SLy4"

# First run: converged HF ground state, written to file for the spectrum
# continuation. The spectrum run needs the file's spenergies as the lower
# spectral bound of the Chebyshev count.
cat << EOF > mocca.data
&nucleus
neutrons=20, protons=20
energy_prec=1e-10
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
nwn = 40, nwp = 40
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='hf'
/
&MomentParam
/
&Cranking
/
EOF

# Second run: frozen-potential spectrum on the same mesh. The spwf set is
# sized by the Chebyshev count to Fermi + 20 MeV (plus one batch of
# overallocation), solved with the fixed-block LOBPCG, straddle-checked and
# trimmed to the cutoff before the output file is written.
cat << EOF > mocca.spectrum.data
&nucleus
neutrons=20, protons=20
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
locking_treshold=1e-6
/
&scfiteration
scfscheme='SPECTRUM'
spectrum_ecut=20.0
/
&wfs
nwn = 40, nwp = 40
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='hf'
OutputFilename='trash'
/
&MomentParam
/
&Cranking
/
EOF

# Run the HF calculation from scratch
./$exe < mocca.data > $outfile.a
mocca_check=$?

# ... and the spectrum run continuing from its wavefunctions
./$exe < mocca.spectrum.data > $outfile.b
mocca_spectrum_check=$?

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Checking
# a) The spectrum run went through the Chebyshev sizing
grep -q 'SPECTRUM SIZING' $outfile.b
check_sizing=$?
# b) It reached the certified straddle check and trimmed to the cutoff
grep -q 'Trimming output spectrum' $outfile.b
check_trim=$?
# c) The trimmed Ritz spectrum reproduces the HF total energy of run 1
Ehf=$(get_total_energy_stdout $outfile.a)
Espectrum=$(get_total_energy_stdout $outfile.b)
compare_floats $Espectrum $Ehf 0.001
check_energy=$?

# remove working directory and traces of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
exit $(($mocca_check || $mocca_spectrum_check || $check_sizing || $check_trim || $check_energy ))
