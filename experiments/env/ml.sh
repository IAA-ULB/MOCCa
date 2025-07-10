#!/bin/bash
# generic module loader for 
# - the Tier-2 machines of CalcUA
# - LUMI

#-------------------------------------------------------------------------------
# FUNCTIONS
#-------------------------------------------------------------------------------

#-------------------------------------------------------------------------------
# CheckSource - Check that this script is sourced. Emit an error message if not.
CheckSource()
{
if [[ $_ == $0 ]]
then
    echo "ERROR: This script must be sourced!" >&2
    exit 1
fi
}

#-------------------------------------------------------------------------------
# on_vaughan - return 1 if run on a Tier-2 cluster of CalcUA, 0 otherwise.
on_vaughan () {
    if [[ "${VSC_INSTITUTE}" == antwerpen ]]; then
        # we are on a Tier-2 cluster of Antwerpen
        return 1
    else
        # we are NOT on a Tier-2 cluster of Antwerpen
        return 0
    fi;
}

#-------------------------------------------------------------------------------
# on_lumi - return 1 if run on a LUMI, 0 otherwise.
on_lumi () {
    # test if environment variable LUMI_MODULEPATH_ROOT has zero length
    if [[ -z "${LUMI_MODULEPATH_ROOT}" ]]; then
        # LUMI_MODULEPATH_ROOT is not defined, so we conclude that we are NOT 
        # on LUMI.
        return 0
    else
        # LUMI_MODULEPATH_ROOT is defined, so we conclude that we are on LUMI.
        return 1
    fi;
}

#-------------------------------------------------------------------------------
# MAIN PROGRAM
#-------------------------------------------------------------------------------
on_vaughan
if [ $? == 1 ]; then
    echo "on_vaughan : yes"
    . ./vaughan/ml.sh -p -v
    return 0
fi

on_lumi
if [ $? == 1 ]; then
    echo "on_lumi : yes"
    . ./lumi/ml.sh
    ml
    return 0
fi

echo "ERROR: Failed to find out where this script was run."
return 1

