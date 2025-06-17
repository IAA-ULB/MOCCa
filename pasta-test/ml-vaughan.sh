#! /bin/bash 
# Use as:
# > . ./ml-vaughan.sh [-v]
# The option -v selects listing currently loaded modules.

module --force purge

# ml calcua/2024a # openmpi seems to be broken 

# the test works fine with these modules:
# ml calcua/2023a
# ml ScaLAPACK

# the test works fine with these modules:
ml calcua/2024a
ml intel
ml iimkl

if [[ "$1" = "-v" ]]; then 
    ml
fi