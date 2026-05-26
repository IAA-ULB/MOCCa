#-------------------------------------------------------------------------------
# Memory leak integration test: verify that memory usage is independent of
# iteration count within a given tolerance.
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# This script tests:
#  - Memory usage (maximum resident size) should not grow significantly with
#    iteration count
#
# Useage
# ------
#   bash memory_leak.sh [EXESUFFIX] [-v|--verbose]
#
# where EXESUFFIX specifies the suffix of the executable to be used
#            Example: "BXL" for "Tantalus.BXL.exe".
#
#   -v, --verbose: Report memory used as a function of iteration count in a table
#
# Dependencies: /usr/bin/time, bc
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner: Memory leak detection test
#-------------------------------------------------------------------------------

set -e

# Parse arguments
VERBOSE=false
EXESUFFIX=""

# Check for verbose flag
for arg in "$@"; do
    case "$arg" in
        -v|--verbose)
            VERBOSE=true
            ;;
        *)
            if [ -z "$EXESUFFIX" ]; then
                EXESUFFIX="$arg"
            else
                echo "Unknown argument: $arg"
                exit 1
            fi
            ;;
    esac
done

# Tolerance for memory comparison (in KB) - memory should not grow by more than this
MEM_TOLERANCE_KB=3000  # 3 MB tolerance by default

# Iteration counts to test
ITERATIONS=(100 500 1000)

# Arrays to store results
declare -a MEMORY_USAGE
declare -a ITERATION_COUNT

export MKL_NUM_THREADS=1

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

echo "========================================="
echo "  Memory Leak Integration Test"
echo "========================================="
echo ""

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Function to run a single test with given maxiter and return memory usage
run_memory_test() {
    local maxiter=$1
    local out_file=$2
    local time_file
    local max_rss

    # Set up test environment
    setup_test_env "memory_leak_${maxiter}" "$EXESUFFIX" "SLy4"

    # Create runtime data with energy_prec=1e-20 to prevent early stopping
    cat <<EOF >tant.data
&nucleus
neutrons=8, protons=8,
energy_prec=1e-20
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
maxiter=${maxiter}
/
&scfiteration
/
&wfs
nwn = 15, nwp = 15
osc_freq = 0.2, 0.2, 0.2
/
&IO
InputFilename='init'
OutputFilename='trash'
/
&MomentParam
/
&Cranking
/
EOF

    # Create temp file for time output
    time_file=$(mktemp)

    # Run the calculation with /usr/bin/time -v to get memory stats
    # /usr/bin/time writes its stats to stderr, so we redirect that to our temp file
    # Program's stdout goes to outfile, program's stderr is discarded (or could be logged)
    /usr/bin/time -v ./$exe <tant.data >"$out_file" 2>"$time_file"
    tantalus_check=$?

    # Extract maximum resident set size (in KB) from time output
    max_rss=$(grep "Maximum resident set size" "$time_file" | awk '{print $6}')
    rm -f "$time_file"

    # Clean up
    teardown_test_env

    echo "$max_rss"
    return $tantalus_check
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run tests for each iteration count
echo "Running memory tests for different iteration counts..."
echo ""

for maxiter in "${ITERATIONS[@]}"; do
    printf "  Testing with maxiter=%4d... " $maxiter
    
    max_rss=$(run_memory_test $maxiter "../logs/memory_leak_${maxiter}.${EXESUFFIX}.out")
    check=$?
    
    if [ $check -ne 0 ]; then
        echo "FAILED (Tantalus exit code: $check)"
        exit 1
    fi
    
    MEMORY_USAGE+=("$max_rss")
    ITERATION_COUNT+=("$maxiter")
    
    echo "Memory: ${max_rss} KB"
done

echo ""

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Verbose output: display table
if [ "$VERBOSE" = true ]; then
    echo "  Memory Usage vs Iteration Count"
    echo "  ----------------------------------------"
    printf "  %-12s | %15s | %10s\n" "Iterations" "Memory (KB)" "Memory (MB)"
    echo "  ----------------------------------------"
    
    for i in ${!ITERATION_COUNT[@]}; do
        iter=${ITERATION_COUNT[$i]}
        mem_kb=${MEMORY_USAGE[$i]}
        mem_mb=$(echo "scale=2; $mem_kb / 1024" | bc)
        printf "  %-12d | %15d | %10.2f\n" $iter $mem_kb $mem_mb
    done
    
    echo "  ----------------------------------------"
    echo ""
fi

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Check that memory usage is independent of iteration count
# Compare memory between first and last test
first_mem=${MEMORY_USAGE[0]}
last_mem=${MEMORY_USAGE[-1]}

mem_diff=$(echo "$last_mem - $first_mem" | bc)
mem_diff_abs=${mem_diff#-}  # Absolute value

echo "Memory difference between ${ITERATIONS[0]} and ${ITERATIONS[-1]} iterations: ${mem_diff_abs} KB"

# Check if difference is within tolerance
if [ 1 -eq "$(echo "$mem_diff_abs <= $MEM_TOLERANCE_KB" | bc)" ]; then
    echo "RESULT: PASSED - Memory usage is stable within tolerance (${MEM_TOLERANCE_KB} KB)"
    exit 0
else
    echo "RESULT: FAILED - Memory grew by ${mem_diff_abs} KB, exceeding tolerance of ${MEM_TOLERANCE_KB} KB"
    exit 1
fi
