#!/bin/bash

# script for building ELPA+pyelpa on lumi
ml LUMI
ml partition/C
ml EasyBuild-user

eb ./ELPA-2024.05.001-cpeGNU-24.03-CPU.eb -r --rebuild