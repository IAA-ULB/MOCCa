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
#   bash memory_leak.sh [EXESUFFIX] [-v|--verbose] [-n|--num-runs NUM_RUNS]
#
# where EXESUFFIX specifies the suffix of the executable to be used
#            Example: "BXL" for "MOCCa.BXL.exe".
#
#   -v, --verbose: Report memory used as a function of iteration count in a table
#   -n, --num-runs: Number of runs per iteration count (default: 1).
#                   The median memory usage across runs is used.
#
# Dependencies: /usr/bin/time, bc, python3
# - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Owner: Memory leak detection test
#-------------------------------------------------------------------------------

set -e

# Parse arguments
VERBOSE=false
EXESUFFIX=""
NUM_RUNS=1

# Check for arguments
while [ $# -gt 0 ]; do
    case "$1" in
        -v|--verbose)
            VERBOSE=true
            shift
            ;;
        -n|--num-runs)
            if [ -n "$2" ] && [ "$2" -eq "$2" ] 2>/dev/null; then
                NUM_RUNS=$2
                shift 2
            else
                echo "Error: -n/--num-runs requires a positive integer argument"
                exit 1
            fi
            ;;
        *)
            if [ -z "$EXESUFFIX" ]; then
                EXESUFFIX="$1"
                shift
            else
                echo "Unknown argument: $1"
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
export OMP_NUM_THREADS=1
#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Basic starting point of all testing scripts
source ../functions.sh

echo "========================================="
echo "  Memory Leak Integration Test"
echo "========================================="
echo ""

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Function to run a single memory test using /usr/bin/time -v
# Returns the Maximum resident set size in KB
run_single_memory_test() {
    local maxiter=$1
    local out_file=$2
    local time_file
    local max_rss

    # Set up test environment
    setup_test_env "memory_leak_${maxiter}" "$EXESUFFIX" "SLy4"

    # Create runtime data with energy_prec=1e-20 to prevent early stopping
    cat <<EOF >mocca.data
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

    # Run with /usr/bin/time -v, sending MOCCa stdout to file and capturing time stats from stderr
    time_file=$(mktemp)
    /usr/bin/time -v ./$exe <mocca.data >"$out_file" 2>"$time_file"
    mocca_check=$?
    
    # Extract Maximum resident set size (in KB) from time output file
    max_rss=$(grep "Maximum resident set size" "$time_file" | awk '{print $6}')
    rm -f "$time_file"

    # Clean up
    teardown_test_env

    if [ -z "$max_rss" ]; then
        echo "ERROR: Could not extract memory usage from time output" >&2
        return 1
    fi

    echo "$max_rss"
    return $mocca_check
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Function to run multiple tests and return the median memory usage
run_memory_test() {
    local maxiter=$1
    local out_file_base=$2
    local -a memories
    local median_mem

    for i in $(seq 1 $NUM_RUNS); do
        local out_file="${out_file_base}.${i}"
        local mem
        mem=$(run_single_memory_test $maxiter "$out_file")
        local check=$?
        
        if [ $check -ne 0 ]; then
            return $check
        fi
        memories+=("$mem")
    done

    # Sort and take median
    IFS=$'\n' sorted=($(printf "%s\n" "${memories[@]}" | sort -n))
    median_mem=${sorted[$((NUM_RUNS / 2))]}
    
    echo "$median_mem"
    return 0
}

#- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
# Run tests for each iteration count
echo "Running memory tests for different iteration counts (${NUM_RUNS} runs each, taking median)..."
echo ""

for maxiter in "${ITERATIONS[@]}"; do
    printf "  Testing with maxiter=%4d... " $maxiter
    
    max_rss=$(run_memory_test $maxiter "../logs/memory_leak_${maxiter}.${EXESUFFIX}")
    check=$?
    
    if [ $check -ne 0 ]; then
        echo "FAILED (MOCCa exit code: $check)"
        exit 1
    fi
    
    MEMORY_USAGE+=("$max_rss")
    ITERATION_COUNT+=("$maxiter")
    
    echo "Memory (median of ${NUM_RUNS}): ${max_rss} KB"
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
# Check that memory usage is independent of iteration count using linear fit
# Perform linear regression: memory = a + b * iterations
# If slope b is close to zero (within tolerance), memory is stable

# Prepare data for Python linear regression
iters_json=$(printf "%s," "${ITERATION_COUNT[@]}")
mems_json=$(printf "%s," "${MEMORY_USAGE[@]}")

# Use Python to compute linear regression slope
slope=$(python3 - <<EOF
import sys
iters = [${iters_json%,}]
mems = [${mems_json%,}]
n = len(iters)
sum_x = sum(iters)
sum_y = sum(mems)
sum_xy = sum(x * y for x, y in zip(iters, mems))
sum_x2 = sum(x ** 2 for x in iters)
denominator = n * sum_x2 - sum_x ** 2
if denominator == 0:
    print("0")
else:
    slope = (n * sum_xy - sum_x * sum_y) / denominator
    print(f"{slope:.6f}")
EOF
)

# Convert slope to integer KB per iteration for comparison
slope_kb_per_iter=$(echo "scale=2; $slope + 0" | bc)

# Check if slope is within tolerance (converted to per-iteration)
# MEM_TOLERANCE_KB is total tolerance over the range, so per-iteration:
# tolerance_per_iter = MEM_TOLERANCE_KB / (max_iter - min_iter)
min_iter=${ITERATIONS[0]}
max_iter=${ITERATIONS[-1]}
iter_range=$(echo "$max_iter - $min_iter" | bc)
tolerance_per_iter=$(echo "scale=6; $MEM_TOLERANCE_KB / $iter_range" | bc)

echo "Linear fit slope: ${slope_kb_per_iter} KB/iteration"
echo "Tolerance: ${tolerance_per_iter} KB/iteration (${MEM_TOLERANCE_KB} KB over ${iter_range} iterations)"

# Check if absolute value of slope is within per-iteration tolerance
slope_abs=$(echo "${slope_kb_per_iter#-}" | bc)
if (( $(echo "$slope_abs <= $tolerance_per_iter" | bc) )); then
    echo "RESULT: PASSED - Memory growth rate (${slope_kb_per_iter} KB/iter) is within tolerance"
    exit 0
else
    echo "RESULT: FAILED - Memory grows at ${slope_kb_per_iter} KB/iteration, exceeding tolerance of ${tolerance_per_iter} KB/iteration"
    exit 1
fi
