#!/bin/bash

on_vaughan () {
    if [[ "${VSC_INSTITUTE}" == antwerpen ]]; then
        return 1
    else
        return 0
    fi;
}

on_lumi () {
    if [[ "${HOME}" == /users/entijske ]]; then
        return 1
    else
        return 0
    fi;
}

on_vaughan
if [ $? == 1 ]; then
    echo on_vaughan
    . ./vaughan/ml.sh -p -v
    return 0
fi

on_lumi
if [ $? == 1 ]; then
    echo on_lumi
    . ./lumi/ml.sh
    ml
    return 0
fi

