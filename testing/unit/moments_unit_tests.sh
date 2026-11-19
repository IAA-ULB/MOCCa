#!/usr/bin/env sh
#-------------------------------------------------------------------------------
# Run the moments unit tests for an O16 nucleus in a minimal box, for three
#  symmetry configurations describing the same physical box: the
#  time-reversal- and parity-conserving build, the time-reversal-broken build
#  and the parity-broken build. The per-run tests are defined in the Fortran
#  code of the moments_testing module: a dense dump of all multipole matrices
#  with l <= 4, the trace check of those matrices against the many-body moment
#  values, and the full-box reconstruction oracle - the reference path for
#  the imaginary operator half (design summary, amended step 4 of section
#  10), which recomputes selected elements by independent mesh quadrature of
#  the chamber-extended spinors and checks identity, STy characteristic and
#  Hermiticity against physics.
#
#  On top of the three runs, the script cross-compares the runs: for every
#  multipole operator that all three builds represent - the even-m operators,
#  (odd m is structurally zero in the TR-conserving build) - the traces of
#  the matrices with the diagonal density,
#      Tr( rho Q ) = sum_i rho_can(i) < i | Q | i >,
#  printed by the three executables must agree: they are the basis-independent
#  carriers of the matrix content, and the three builds store the same
#  operator in different symmetry-adapted bases of the same state.
#  Element-by-element equality of the matrices is meaningless (different
#  bases); the traces are the invariants that must not differ beyond build
#  differences. The tests are NOT a convergence test: the runs are started
#  from the same initialisation (which is identical for the three builds by
#  construction) and take zero SCF iterations, so the compared state is
#  independent of the iteration path.
#
#  ATTENTION on the input files: the physical box is 24 fm in every direction.
#  The z-axis of the conservative and TR-broken runs is folded (dv carries
#  the multiplicity), so nz=12 covers the full length there; the parity-broken
#  run stores the z-axis unfolded and hence needs nz=24 to describe the same
#  box. The same O16/SLy4 nucleus should converge to the same many-body state
#  in all three runs - a large discrepancy in the final total energies is a
#  sign that the inputs no longer describe the same box or state.
#
#  The unit tests are armed with unit_test=.true. in the &mesh/ namelist: each
#  run executes the tests on the final mean-field state and stops with the
#  exit status of the test suite.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Useage
# ------
#   bash moments_unit_tests.sh [SUFFIX [SUFFIX_T SUFFIX_P]]
#         [--pairing HF|BCS|HFB] [--maxiter N]
#
# where
# SUFFIX           : executable suffix of the conservative build, e.g. "NLO"
#                    for MOCCa.NLO.exe (default: NLO)
# SUFFIX_T, SUFFIX_P: the TR-broken and parity-broken executables. Default:
#                    derived from SUFFIX (MOCCa.NLO-T.exe / MOCCa.NLO-P.exe);
#                    CI passes them explicitly, since its suffix scheme puts
#                    the configuration first (MOCCa.NLO-T.gfortran.NUCLEI...).
# --pairing       : pairing type of the mean-field state (default: HF).
#                   Note that the trace check is skipped for HFB, whose
#                   density matrix is canonical.
# --maxiter       : number of SCF iterations per run (default: 0 - the runs
#                   act on the common initialisation; see above). Kept for
#                   debugging runs that iterate towards self-consistency.
#
# Dependencies: none - but the three executables must exist.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner                : W. Ryssens [wouter.ryssens@ulb.be]
# Reference commit hash: this script lands with the moments_testing module
#-------------------------------------------------------------------------------
# Initialize default values
suffix="NLO"     # Default base configuration
pairing="HF"     # Default pairing type
maxiter=0        # SCF iterations per run: zero - this is not a convergence
                 #  test, the identical initialisation is the whole point

# Tolerances of the cross-build comparison: the traces of the common
#  operators must agree within  atol + rtol*max|value|  across the three runs.
#  The runs agree to all printed digits when their symmetry-adapted boxes
#  coincide (NLO vs NLO-T) and to ~1e-8 absolute for the unfolded-z
#  parity-broken run; a real machinery error changes these traces at the
#  O(1) relative level.
crtol=1.0e-6
catol=1.0e-7

# Temporary array to hold arguments
args=()

while [[ $# -gt 0 ]]; do
  case "$1" in
    --pairing)
      shift
      if [[ "$1" == "HF" || "$1" == "BCS" || "$1" == "HFB" ]]; then
        pairing="$1"
      else
        echo "Error: --pairing must be HF, BCS or HFB"
        exit 1
      fi
      shift
      ;;
    --maxiter)
      shift
      if [[ "$1" =~ ^[0-9]+$ ]]; then
        maxiter="$1"
      else
        echo "Error: --maxiter must be a positive integer"
        exit 1
      fi
      shift
      ;;
    *)
      args+=("$1")
      shift
      ;;
  esac
done

# Restore the positional arguments: the executable suffix of the conservative
#  build, optionally followed by those of the TR-broken and parity-broken
#  builds. By default the latter two are derived from the conservative suffix
#  (the local convention MOCCa.NLO-T.exe); CI hands them over explicitly.
set -- "${args[@]}"
if [[ "$1" != "" ]]; then
  suffix="$1"
fi
suffixT="${2:-$suffix-T}"
suffixP="${3:-$suffix-P}"

#- - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

# Set up: copies MOCCa.$suffix.exe + the parameterization; the broken-symmetry
#  executables are added by hand.
# ATTENTION: this script lives in testing/unit, one level deeper than the
#  testing/integration scripts that setup_test_env's relative EXECDIR is
#  written for: correct it for use after the cd into work/.
setup_test_env "moments_unit_tests" "$suffix" "SLy4"
EXECDIR="../../../exec/"
if [ ! -f "$EXECDIR/MOCCa.$suffixT.exe" ] || [ ! -f "$EXECDIR/MOCCa.$suffixP.exe" ]; then
  echo "Error: this test needs MOCCa.$suffix.exe, MOCCa.$suffixT.exe and MOCCa.$suffixP.exe."
  exit 1
fi
cp $EXECDIR/MOCCa.$suffixT.exe .
cp $EXECDIR/MOCCa.$suffixP.exe .
exeT="MOCCa.$suffixT.exe"
exeP="MOCCa.$suffixP.exe"
outfileT="../logs/moments_unit_tests.$suffixT.out"
outfileP="../logs/moments_unit_tests.$suffixP.out"

# - - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - -
# Input files: the same O16/SLy4 nucleus in the same physical 24 fm box for
#  all three runs. The TR-broken and parity-broken builds demand even numbers
#  of spwfs and allowtransform=.true. in &IO/ (the goodqn/transform recipes);
#  the same even nwn/nwp is used for the conservative run for comparability.
#  The parity-broken run unfolds the z-axis: nz=24 covers the box that the
#  folded runs cover with nz=12.
for run in "mocca.data" "moccaT.data" "moccaP.data"; do
  extra=""
  nzval=12
  if [ "$run" = "moccaT.data" ]; then
    extra="allowtransform=.true."
  fi
  if [ "$run" = "moccaP.data" ]; then
    extra="allowtransform=.true."
    nzval=24
  fi
  cat << EOF > $run
&nucleus
neutrons=8, protons=8
/
&mesh
nx=12, ny=12, nz=$nzval, dx=1.0
unit_test=.true.
/
&func
name_param='SLy4'
/
&pairing
type="$pairing"
/
&evolution
maxiter=$maxiter
/
&scfiteration
/
&wfs
nwn = 16, nwp = 16
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
$extra
/
&MomentParam
/
&Cranking
/
EOF
done

# - - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - -
# The three runs. Each exits with the status of its own test suite.
./$exe   < mocca.data   > $outfile
mocca_check=$?
./$exeT  < moccaT.data  > $outfileT
moccaT_check=$?
./$exeP  < moccaP.data  > $outfileP
moccaP_check=$?

# - - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - -
# Cross-build comparison of the common operators: parse the trace tables of
#  the three logs (fixed-column rows "l m half ...", scoped by the section
#  banner), keep the rows printed in all three runs with even m, and compare
#  the trace columns. The final total energies are printed for context: a
#  large discrepancy there means the runs no longer describe the same state
#  and the comparison is void.
# ATTENTION: LC_ALL=C is mandatory for the number parsing (locale warning).
cross_check=0
LC_ALL=C awk -v rtol=$crtol -v atol=$catol '
FNR==1 { sec=""; f++ }
/Trace of the multipole matrices/            { sec="trace" ; next }
/largest deviation/                          { sec=""      ; next }
/Total energy:/                              { energy[f] = $3 }
sec=="trace" && $1 ~ /^[0-9]+$/ && $2 ~ /^[0-9]+$/ && ($3=="Re" || $3=="Im") {
  key = $1 "_" $2 "_" $3
  tr[key "," f ",n"] = $4
  tr[key "," f ",p"] = $6
}
END {
  printf("final total energies: %.6f / %.6f / %.6f MeV\n", energy[1], energy[2], energy[3])
  nrow = 0 ; failed = 0 ; worst = 0
  for (keyf in tr) {
    split(keyf, a, ",")
    if (a[2] != 1 || a[3] != "n") continue             # each key handled once
    key = a[1]
    if ( !((key ",2,n") in tr) || !((key ",3,n") in tr) ) continue
    split(key, q, "_")
    if (q[2] % 2 != 0) continue                        # common = even m only
    nrow++
    for (isonum=1; isonum<=2; isonum++) {
      iso = (isonum==1) ? "n" : "p"
      x1 = tr[key ",1," iso] ; x2 = tr[key ",2," iso] ; x3 = tr[key ",3," iso]
      mx = (x1>x2 ? (x1>x3?x1:x3) : (x2>x3?x2:x3))
      dev = (x1>x2 ? x1-x2 : x2-x1)
      d13 = (x1>x3 ? x1-x3 : x3-x1) ; if (d13 > dev) dev = d13
      d23 = (x2>x3 ? x2-x3 : x3-x2) ; if (d23 > dev) dev = d23
      tol = atol + rtol*(mx<0 ? -mx : mx)
      if (dev/tol > worst) worst = dev/tol
      if (dev > tol) { failed=1 ; printf("  trace row %s (%s): %.10E %.10E %.10E\n", key, iso, x1, x2, x3) }
    }
  }
  printf("cross-build comparison: %d common operator rows, worst tolerance usage %.2f\n", nrow, worst)
  exit failed
}' $outfile $outfileT $outfileP
cross_check=$?

teardown_test_env

#- - - - - - - - - - - - - - - -  -- - - - - - - - - - - - - - - - - - - - - -
echo ""
echo "========= Test Summary ========"
echo "NLO  moments unit tests status : $mocca_check"
echo "NLO-T moments unit tests status: $moccaT_check"
echo "NLO-P moments unit tests status: $moccaP_check"
echo "cross-build comparison status  : $cross_check"
echo "logs at testing/unit/logs/moments_unit_tests.$suffix.out , $suffixT.out , $suffixP.out"
echo "==============================="

[ $mocca_check -eq 0 ] && [ $moccaT_check -eq 0 ] && [ $moccaP_check -eq 0 ] && [ $cross_check -eq 0 ]
exit $?
