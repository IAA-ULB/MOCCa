#-------------------------------------------------------------------------------
# Debug dumpfile test: run a minimal O16 calculation with writedumpfile=.true.
# and verify the resulting <basename>.dump.hdf5 file:
#  - it exists and contains the expected datasets;
#  - <psi|h|psi> as computed from the stored hpsi matches the spenergies
#    printed in the .hdf5 output file for a few spwfs.
#
# Useage
# ------
#   bash dumpfile.sh [EXESUFFIX]
#
# where EXESUFFIX specifies the suffix of the executable to be used
#            Example: "NLO" for "MOCCa.NLO.exe".
#
# Dependencies: minimal.sh input; python3 with h5py for the dataset checks.
#-------------------------------------------------------------------------------

set -e
#- - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up
setup_test_env "dumpfile" "$1" "SLy4"

# Create runtime data: the minimal O16 calculation, with the dumpfile switched on
cat << EOF > mocca.data
&nucleus
neutrons=8, protons=8
writedumpfile=.true.
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
OutputFilename='mf.hdf5'
/
&MomentParam
/
&Cranking
/
EOF

# Run the calculation
./$exe < mocca.data > $outfile
mocca_check=$?

#- - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - -
# a) The dumpfile exists and contains the expected datasets
dumpcheck=$(python3 ../dumpfile_check.py mf.dump.hdf5 mf.hdf5 2>&1)
dump_status=$?
echo "$dumpcheck"
if [ $dump_status -ne 0 ]; then
  echo "Dumpfile check FAILED."
  teardown_test_env
  exit 1
fi

# remove working directory and traces of these calculations
teardown_test_env
#- - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - - -
# Return exit code 1 if any of the checks failed
exit $mocca_check
