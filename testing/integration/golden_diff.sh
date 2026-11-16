#! /usr/bin/env bash
#-------------------------------------------------------------------------------
# Golden-output regression: rerun an integration test and compare its STDOUT,
# line by line and number by number, against a golden reference output that
# was produced by a trusted (i.e. un-refactored) executable.
#
# How golden references are made
# ------------------------------
# Run the test you want a golden for with the reference executable. Following
# the setup_test_env convention (see ../functions.sh), a reference for
# 'minimal', based on the NLO reference executable, is produced by
#
#       bash minimal.sh NLO.reference
#       -> writes logs/minimal.NLO.reference.out
#
# i.e. golden references are ordinary test logs with suffix
# '<SUFFIX>.reference'. They are kept under version control: .gitignore
# excludes testing/*/logs/* except *.reference.out.
#
# What is compared
# ----------------
# Everything in the STDOUT except two volatile blocks, which are stripped at
# comparison time (the stored reference itself is always raw):
#   - the version/commit banner (changes with every commit)
#   - the 'Timers of rank...' table at the end (changes every run)
# Every number on the remaining lines must satisfy
#       | a - b |  <=  atol + rtol * max(|a|, |b|)
# so that last-digit drift from reshuffled summation order does not fail the
# comparison, while real structural changes (e.g. sign errors in phase tables,
# altered selection rules) do. Non-numeric tokens, the number of lines and the
# number of tokens per line must match exactly.
#
# Usage
# ------
#   bash golden_diff.sh [TEST] [WORKSUFFIX] [REFSUFFIX] [flags]
#
#   TEST        integration test to run (X.sh)         [default: minimal]
#   WORKSUFFIX  executable under test (MOCCa.X.exe)    [default: NLO]
#   REFSUFFIX   golden reference suffix                [default: <WORKSUFFIX>.reference]
#
#   flags:
#     --no-run    compare logs/<TEST>.<WORKSUFFIX>.out as it is (no rerun)
#     --update    (re)make the current output the golden reference
#     --rtol F    relative tolerance                   [default: 1e-8]
#     --atol F    absolute tolerance                    [default: 1e-7]
#
#   Typical use during the refactor:
#       bash golden_diff.sh minimal NLO
#   which reruns minimal.sh NLO and then compares logs/minimal.NLO.out
#   against logs/minimal.NLO.reference.out.
#
# Notes
# -----
# - The integral checks performed by <TEST>.sh itself (e.g. the 1 keV energy
#   anchor inside minimal.sh) are honoured through its exit code: this script
#   fails when they fail, even if the golden comparison passes.
# - Line numbers in mismatch reports refer to the normalized output (version
#   banner and timers table stripped), not to the raw file.
# - Currently compares the mean-field STDOUT only, i.e. files named
#   logs/<TEST>.<SUFFIX>.out. Tests that also produce .fam.out logs need a
#   small extension of the comparison below.
#
# Dependencies: <TEST>.sh and the ../functions.sh conventions.
#-------------------------------------------------------------------------------

#- - - - - - - - - - -  -- Argument parsing - - - - - - - - - - - - - - - - -
testname=minimal
suffix=NLO
refsuffix=""
norun=0
update=0
rtol=1e-8
atol=1e-7

positional=0
while [ $# -gt 0 ]; do
  case "$1" in
    --no-run)  norun=1 ;;
    --update)  update=1 ;;
    --rtol)    rtol="$2"; shift ;;
    --atol)    atol="$2"; shift ;;
    -*)        echo "golden_diff.sh: unknown flag '$1'" >&2 ; exit 2 ;;
    *)         positional=$((positional+1))
               case $positional in
                 1) testname="$1" ;;
                 2) suffix="$1" ;;
                 3) refsuffix="$1" ;;
               esac ;;
  esac
  shift
done
if [ -z "$refsuffix" ]; then
  refsuffix="${suffix}.reference"
fi

if [ ! -f "../functions.sh" ]; then
  echo "golden_diff.sh: run from inside testing/integration/, like the other tests" >&2
  exit 2
fi

testscript="${testname}.sh"
current="logs/${testname}.${suffix}.out"
reference="logs/${testname}.${refsuffix}.out"

echo "==> golden regression for ${testname} (current: ${suffix}, reference: ${refsuffix})"
echo "==> current  : ${current}"
echo "==> reference: ${reference}"

#- - - - - - - - - - -  -- Run the test - - - - - - - - - - - - - - - - - - -
testrc=0
if [ $norun -eq 0 ]; then
  if [ ! -f "$testscript" ]; then
    echo "golden_diff.sh: no ${testscript} in this directory" >&2
    exit 2
  fi
  echo "==> running: bash ${testscript} ${suffix}"
  bash "$testscript" "$suffix"
  testrc=$?
  echo "==> ${testscript} exited ${testrc}"
  if [ $testrc -ne 0 ] && [ -d work ]; then
    echo "    note: a work/ directory was left behind, inspect or remove it"
  fi
else
  echo "==> --no-run: comparing existing output"
fi

#- - - - - - - -  -- (Re)make the golden reference - - - - - - - - - - - - - -
if [ $update -eq 1 ]; then
  if [ ! -f "$current" ]; then
    echo "golden_diff.sh: nothing to promote, ${current} does not exist" >&2
    exit 2
  fi
  cp "$current" "$reference"
  echo "==> ${reference} refreshed from ${current}"
  echo "    commit the new golden deliberately."
  exit 0
fi

#- - - - - - - - - - -  -- Sanity of all inputs - - - - - - - - - - - - - - -
if [ ! -f "$reference" ]; then
  echo "golden_diff.sh: missing golden reference ${reference}" >&2
  echo "  produce one with the trusted executable, e.g.:" >&2
  echo "      bash ${testscript} ${refsuffix}" >&2
  exit 2
fi
if [ ! -f "$current" ]; then
  echo "golden_diff.sh: ${current} not found - did '${testscript} ${suffix}' run?" >&2
  exit 2
fi

#- - - - - -  -- Normalisation and numeric comparison - - - - - - - - - - - -
# Strip the two volatile blocks: everything between the version banner and
# the symmetry banner, and everything from the timers table to the end.
normalize_stdout() {
  LC_ALL=C awk '
    /Timers of rank/                { exit }
    /Version Information/           { skip = 1 }
    skip && /Symmetry Information/  { skip = 0 }
    !skip { print }
  ' "$1"
}

tmpdir=$(mktemp -d)
trap 'rm -rf "$tmpdir"' EXIT
normalize_stdout "$reference" > "$tmpdir/ref.norm"
normalize_stdout "$current"   > "$tmpdir/cur.norm"

# LC_ALL=C: numbers are parsed and printed with '.' as the decimal separator,
# independent of the user's locale (cf. the LC_NUMERIC guard in get_strength,
# ../functions.sh).
compare_rc=0
LC_ALL=C awk -v ref="$tmpdir/ref.norm" -v cur="$tmpdir/cur.norm" \
    -v rtol="$rtol" -v atol="$atol" '
  # Numeric value of a token, or "" if it is not a number. Punctuation at
  # the edges of a token is ignored; Fortran D exponents are supported.
  function num(tok,    s) {
    s = tok
    sub(/^[^0-9A-Za-z+-]+/, "", s)
    sub(/[^0-9A-Za-z+-]+$/, "", s)
    gsub(/[Dd]/, "E", s)
    if (s ~ /^[+-]?([0-9]+\.?[0-9]*|\.[0-9]+)([Ee][+-]?[0-9]+)?$/)
      return s + 0
    return ""
  }
  BEGIN {
    rn = 0 ; while ((getline l < ref) > 0) R[++rn] = l ; close(ref)
    cn = 0 ; while ((getline l < cur) > 0) C[++cn] = l ; close(cur)

    shown = 0        # mismatches printed so far
    nnum = 0         # numbers compared
    worst = 0 ; worstline = 0 ; worstref = "" ; worstcur = ""
    bad = 0

    common = (rn < cn) ? rn : cn
    for (i = 1 ; i <= common ; i++) {
      lineR = R[i] ; sub(/^[ \t]+/, "", lineR) ; sub(/[ \t]+$/, "", lineR)
      lineC = C[i] ; sub(/^[ \t]+/, "", lineC) ; sub(/[ \t]+$/, "", lineC)
      nR = split(lineR, TR, /[ \t]+/)
      nC = split(lineC, TC, /[ \t]+/)

      if (nR != nC) {
        bad++
        if (shown < 5) {
          shown++
          printf "  MISMATCH line %d: %d tokens (ref) vs %d tokens (cur)\n", i, nR, nC
          printf "    ref: %s\n    cur: %s\n", R[i], C[i]
        }
      }
      for (j = 1 ; j <= nR && j <= nC ; j++) {
        a = num(TR[j]) ; b = num(TC[j])
        if (a != "" && b != "") {
          nnum++
          d = a - b ; if (d < 0) d = -d
          m = (a < 0) ? -a : a ; mb = (b < 0) ? -b : b ; if (mb > m) m = mb
          t = atol + rtol * m
          if (d / t > worst) {
            worst = d / t ; worstline = i ; worstref = TR[j] ; worstcur = TC[j]
          }
          if (d > t) {
            bad++
            if (shown < 5) {
              shown++
              printf "  MISMATCH line %d, number %d: %s (ref) vs %s (cur)\n", i, j, TR[j], TC[j]
              printf "    ref: %s\n    cur: %s\n", R[i], C[i]
            }
          }
        } else if (TR[j] != TC[j]) {
          bad++
          if (shown < 5) {
            shown++
            printf "  MISMATCH line %d, token %d: %s (ref) vs %s (cur)\n", i, j, TR[j], TC[j]
            printf "    ref: %s\n    cur: %s\n", R[i], C[i]
          }
        }
      }
    }
    if (rn != cn) {
      bad++
      printf "  MISMATCH line count: %d (ref) vs %d (cur)\n", rn, cn
      if (rn > cn) printf "    first extra reference line: %s\n", R[cn+1]
      if (cn > rn) printf "    first extra current line:   %s\n", C[rn+1]
    }

    printf "  lines compared: %d   numbers compared: %d\n", common, nnum
    if (worst == 0)
      printf "  all numbers identical\n"
    else
      printf "  worst number used %.3g%% of the tolerance (line %d: %s vs %s)\n", 100*worst, worstline, worstref, worstcur
    exit (bad > 0)
  }
'
compare_rc=$?

#- - - - - - - - - - -  -- Verdict - - - - - - - - - - - - - - - - - - - - -
if [ $compare_rc -eq 0 ] && [ $testrc -eq 0 ]; then
  echo "GOLDEN PASS: ${testname} (${current} vs ${reference})"
  exit 0
else
  echo "GOLDEN FAIL: ${testname} (${current} vs ${reference})"
  [ $testrc -ne 0 ] && echo "  also: ${testscript} itself exited ${testrc}"
  exit 1
fi
