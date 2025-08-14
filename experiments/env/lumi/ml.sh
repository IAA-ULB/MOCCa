#! /bin/bash 

module --force purge

# the test works fine with these modules:
ml LUMI/24.03
ml partition/C

ml PrgEnv-cray
ml buildtools
ml cray-python
ml EasyBuild-user

ml cray-python/3.11.7

alias wq='watch -n 5 squeue -u entijske'