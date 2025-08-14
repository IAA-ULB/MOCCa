 This file must be sourced

module purge

# default modules
ml LUMI
ml partition/C
ml PrgEnv-cray
ml buildtools
ml cray-python
ml cray-libsci


for cla in "$@"
do
  if [[ $cla -eq -p ]]
  then
    # load modules for performance analysis
    ml perftools-base
    ml perftools-lite
  fi
done

ml